#!/usr/bin/env python
"""低丰度转录本 smFISH 探针设计的批量测试（第二轮）。

与第一轮（batch_smfish_test.py）的区别——实验目的变了：
    第一轮  物种区分探针 → 用另外两个物种做背景过滤，考察交叉反应淘汰；
    本轮    每种菌独立液培、独立检测**低丰度转录本**，不存在跨种背景，
            因此**关闭跨物种背景过滤**（hosts = 空列表）。

低丰度检测的参数取向：
    - 目标探针数从 20 提高到 48：低丰度 RNA 的 FISH 信号 ∝ 结合在单个
      转录本上的探针数，可用探针越多检出越可靠（经典 smFISH 对低丰度
      转录本常用 ~48 条探针/转录本）；
    - min_gap = 2 nt：轻微间隔避免相邻探针空间位阻，同时保留堆叠密度；
    - 保留靶标自身基因组的重复序列检测（max_target_hits = 10）：
      这是探针质量的内在 QC，与跨种过滤无关。

"低丰度基因"的选取（注释数据不含表达量，用功能类别做低丰度代理；
严格丰度排序需 RNA-seq 确认）：
    a) σ 因子（sigma factor）
    b) 丝氨酸/苏氨酸蛋白激酶（serine/threonine protein kinase，pkn）
    c) 双组分系统的反应调节因子（response regulator）
    d) 转录调控因子（transcriptional regulator；MycoBrowser 的
       "regulatory proteins" 功能类别）
    以上四类是细菌中经典的低拷贝转录本。每物种按 a→d 顺序取满
    16 个（类内按 locus 排序），三个物种共 ~48 个。

来源：MTB H37Rv 与 M. smegmatis MC2-155 用 MycoBrowser GFF；
BCG Pasteur 1173P2 用 NCBI RefSeq GFF（MycoBrowser 不收录 BCG，
其蛋白编码基因与 H37Rv 一一对应，注释内容等价）。

运行（项目目录）：
    /opt/anaconda3/envs/Probe/bin/python batch_lowabundance_test.py
"""

from __future__ import annotations

import csv
import re
import statistics
import sys
import time
from pathlib import Path

PROJECT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT / "src"))

from batch_smfish_test import (  # noqa: E402  复用第一轮的解析与统计工具
    ANNOT_DIR,
    ORGANISMS,
    extract_gene,
    failure_category,
    load_genome,
    parse_attributes,
    write_fasta,
)
from probedesign.models import DesignParams, DesignResult, ReferenceGenome  # noqa: E402
from probedesign.pipeline import run_design  # noqa: E402

OUT_DIR = PROJECT / "test_data" / "lowabundance_results"
GENE_FASTA_DIR = OUT_DIR / "gene_fastas"
MAX_PER_ORGANISM = 16
DESIRED_PROBES = 48  # 低丰度转录本：尽量多的探针堆叠提高检出

LOW_ABUNDANCE_CLASSES = [
    ("sigma", re.compile(r"sigma factor", re.I)),
    ("kinase", re.compile(r"serine/threonine protein kinase", re.I)),
    ("response_regulator", re.compile(r"response regulator", re.I)),
    ("regulator", re.compile(r"transcriptional regulator", re.I)),
]


def _class_of(gene: str, product: str) -> str | None:
    """返回基因所属的低丰度类别名；不属于任何类别返回 None。"""
    gene_clean = gene.lower()
    if re.fullmatch(r"sig[a-z]", gene_clean):
        return "sigma"
    text = product
    for name, pattern in LOW_ABUNDANCE_CLASSES:
        if pattern.search(text):
            return name
    if re.fullmatch(r"pkn[a-z]", gene_clean):
        return "kinase"
    return None


def select_genes_from_gff(gff_path: Path, *, mycobrowser: bool) -> list[dict]:
    """解析 GFF，按低丰度类别选取基因（类内按 locus 排序，填满 16 个/种）。"""
    buckets: dict[str, list[dict]] = {name: [] for name, _ in LOW_ABUNDANCE_CLASSES}
    with gff_path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attrs = parse_attributes(parts[8])
            if mycobrowser:
                locus = attrs.get("Locus", "")
                gene = attrs.get("Name", "") or locus
                product = attrs.get("Product", "")
                if attrs.get("Is_Pseudogene", "No") == "Yes":
                    continue
            else:
                locus = attrs.get("locus_tag", "")
                gene = attrs.get("gene", "") or locus
                product = attrs.get("product", "")
            cls = _class_of(gene, product)
            if cls is None:
                continue
            buckets[cls].append(
                {
                    "locus": locus,
                    "gene": gene,
                    "product": product[:90],
                    "cls": cls,
                    "start": int(parts[3]),
                    "stop": int(parts[4]),
                    "strand": parts[6],
                }
            )
    # 四类轮转取样（各类内按 locus 排序），让每物种的基因集覆盖
    # sigma / 激酶 / 反应调节因子 / 转录调控因子 四类低丰度转录本。
    buckets = {name: sorted(bucket, key=lambda g: g["locus"])
               for name, bucket in buckets.items()}
    selected: list[dict] = []
    class_cycle = [name for name, _ in LOW_ABUNDANCE_CLASSES]
    index = 0
    while len(selected) < MAX_PER_ORGANISM and any(buckets.values()):
        cls = class_cycle[index % len(class_cycle)]
        index += 1
        if buckets[cls]:
            selected.append(buckets[cls].pop(0))
    return selected


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    GENE_FASTA_DIR.mkdir(parents=True, exist_ok=True)

    # ---- 第 1 步：选基因 ----
    selected = {
        "mtb": select_genes_from_gff(
            ANNOT_DIR / "mycobrowser_h37rv.gff", mycobrowser=True
        ),
        "msm": select_genes_from_gff(
            ANNOT_DIR / "mycobrowser_smegmatis.gff", mycobrowser=True
        ),
        "bcg": select_genes_from_gff(
            ANNOT_DIR / "GCF_000009445.1_ASM944v1_genomic.gff", mycobrowser=False
        ),
    }
    total = sum(len(v) for v in selected.values())
    print("[1/3] 低丰度基因选取完成：", flush=True)
    for key, genes in selected.items():
        classes = [g["cls"] for g in genes]
        counts = {c: classes.count(c) for c in dict.fromkeys(classes)}
        print(f"      {ORGANISMS[key]['label']}: {len(genes)} 个 {counts}", flush=True)
    print(f"      合计 {total} 个", flush=True)

    # ---- 第 2 步：提取序列 ----
    genome_seqs = {key: load_genome(cfg["genome"]) for key, cfg in ORGANISMS.items()}
    gene_meta: list[dict] = []
    for org_key, genes in selected.items():
        for gene in genes:
            seq = extract_gene(genome_seqs[org_key], gene)
            fasta_id = f"{org_key}|{gene['locus']}|{gene['gene']}"
            path = GENE_FASTA_DIR / f"{org_key}_{gene['locus']}.fasta"
            write_fasta(path, fasta_id, seq)
            gene_meta.append({**gene, "org": org_key, "length_nt": len(seq),
                              "fasta_path": str(path), "fasta_id": fasta_id})
    print(f"[2/3] 序列提取完成（{len(gene_meta)} 条），开始批量设计…", flush=True)

    # ---- 第 3 步：批量设计（无跨物种背景；只做自身基因组重复检测）----
    params = DesignParams(
        design_scheme="smFISH",
        desired_probe_count=DESIRED_PROBES,
        min_gap=2,
    )
    rows: list[dict] = []
    results_by_gene: dict[str, DesignResult] = {}
    started = time.time()
    for i, meta in enumerate(gene_meta, 1):
        t0 = time.time()
        try:
            result: DesignResult = run_design(
                meta["fasta_path"], ORGANISMS[meta["org"]]["index"], [], params, threads=4
            )
        except Exception as exc:
            rows.append({**meta, "error": str(exc)[:200]})
            print(f"  ({i}/{len(gene_meta)}) {meta['org']}:{meta['gene']} 失败：{exc}", flush=True)
            continue
        results_by_gene[meta["fasta_id"]] = result

        passed = result.passed_probes
        tms = [p.tm for p in passed]
        thermo_pass = sum(
            1 for p in result.probes if failure_category(p) != "thermo"
        )
        spec_pass = sum(
            1 for p in result.probes if failure_category(p) not in ("thermo", "specificity")
        )
        target_rejected = sum(
            1 for p in result.probes
            if any("target_hits" in r for r in p.failure_reasons)
        )
        rows.append(
            {
                **meta,
                "candidates": len(result.probes),
                "thermo_pass": thermo_pass,
                "spec_pass": spec_pass,
                "final_probes": len(passed),
                "coverage_pct": round(
                    sum(p.stop - p.start for p in passed) / result.target_length * 100, 1
                ) if result.target_length else 0.0,
                "tm_mean": round(statistics.mean(tms), 1) if tms else "",
                "tm_sd": round(statistics.stdev(tms), 1) if len(tms) > 1 else 0.0,
                "target_rejected": target_rejected,
                "runtime_s": round(time.time() - t0, 1),
                "error": "",
            }
        )
        print(
            f"  ({i}/{len(gene_meta)}) {meta['org']}:{meta['gene']}({meta['cls']}) "
            f"候选{rows[-1]['candidates']} → 最终{len(passed)} "
            f"覆盖{rows[-1]['coverage_pct']}% ({rows[-1]['runtime_s']}s)",
            flush=True,
        )

    # ---- 第 4 步：落盘 ----
    with (OUT_DIR / "batch_summary.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    detail_rows: list[dict] = []
    for meta, row in zip(gene_meta, rows):
        if row.get("error") or meta["fasta_id"] not in results_by_gene:
            continue
        result = results_by_gene[meta["fasta_id"]]
        for probe in result.passed_probes:
            detail_rows.append(
                {
                    "org": meta["org"],
                    "locus": meta["locus"],
                    "gene": meta["gene"],
                    "cls": meta["cls"],
                    "probe_id": probe.probe_id,
                    "start": probe.start,
                    "stop": probe.stop,
                    "sequence": probe.sequence,
                    "tm": round(probe.tm, 1),
                    "gc": round(probe.gc_content, 2),
                    "score": round(probe.score, 3),
                }
            )
    with (OUT_DIR / "all_final_probes.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(detail_rows[0].keys()))
        writer.writeheader()
        writer.writerows(detail_rows)

    elapsed = time.time() - started
    total_probes = sum(int(r["final_probes"]) for r in rows if not r.get("error"))
    print(f"[3/3] 完成：{len(gene_meta)} 个基因 / {total_probes} 条探针，"
          f"总耗时 {elapsed / 60:.1f} 分钟", flush=True)
    print(f"      汇总：{OUT_DIR / 'batch_summary.csv'}", flush=True)


if __name__ == "__main__":
    main()
