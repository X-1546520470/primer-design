#!/usr/bin/env python
"""Propionyl-CoA 代谢基因的 smFISH 探针设计批量测试。

目的：
    用真实的分枝杆菌基因集，对 ProbeStudio 的 smFISH 设计工作流做端到端
    批量测试（候选枚举 → 热力学过滤 → bowtie2 靶标/背景比对 → 打分 → 选点）。

基因选取（共 ≥30 个，三个物种）：
    来源 1  MycoBrowser GFF（M. tuberculosis H37Rv、M. smegmatis MC2-155）
    来源 2  NCBI RefSeq GFF（M. bovis BCG Pasteur 1173P2——MycoBrowser 不
            收录 BCG，故按 H37Rv 的 MycoBrowser 功能注释取 BCG 同源基因；
            两基因组同源度 >99.9%，基因一一对应）
    筛选标准（满足其一）：
        a) 功能注释含 propionyl-CoA / methylcitrate / methylisocitrate /
           methylmalonyl-CoA（丙酰辅酶A直接代谢酶）
        b) 基因名属于 prpB/prpC/prpD/prpR（甲基柠檬酸循环）
        c) 基因名属于 icl1/icl2/aceAa/aceAb/aceA（乙醛酸支路——MTB 中
           丙酰辅酶A脱毒的替代通路，与甲基柠檬酸循环功能冗余）
    每个物种最多取 16 个，按 locus 标签排序后取前若干个。

背景基因组设置（交叉反应过滤的方向性）：
    MTB 基因  → 背景 = BCG + MSM
    BCG 基因  → 背景 = MTB + MSM
    MSM 基因  → 背景 = MTB + BCG
    注意：MTB 与 BCG 同源度极高，保守基因的探针会被背景过滤大量淘汰——
    这是工具的预期行为（说明这些基因不适合做区分二者的探针），
    报告中将单独统计该现象。

运行（在项目目录）：
    /opt/anaconda3/envs/Probe/bin/python batch_smfish_test.py
"""

from __future__ import annotations

import csv
import json
import re
import sys
import time
from dataclasses import asdict
from pathlib import Path

# 让脚本在未安装包的情况下也能从 src/ 导入引擎
PROJECT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT / "src"))

from probedesign.models import DesignParams, DesignResult, ReferenceGenome  # noqa: E402
from probedesign.pipeline import run_design  # noqa: E402
from probedesign.utils import reverse_complement  # noqa: E402

ANNOT_DIR = PROJECT / "test_data" / "annotations"
GENOME_DIR = PROJECT / "genomes"
INDICES_DIR = PROJECT / "indices"
OUT_DIR = PROJECT / "test_data" / "batch_results"
GENE_FASTA_DIR = OUT_DIR / "gene_fastas"

# 关键词与基因名的筛选规则（见模块 docstring）
KEYWORDS = re.compile(r"propionyl|methylcitrate|methylisocitrate|methylmalonyl", re.I)
CORE_GENES = {"prpB", "prpC", "prpD", "prpR", "icl1", "icl2", "aceAa", "aceAb", "aceA"}
MAX_PER_ORGANISM = 16

# 三个物种的基因组与索引配置
ORGANISMS = {
    "mtb": {
        "label": "M. tuberculosis H37Rv",
        "genome": GENOME_DIR / "mtb_h37rv.fna",
        "index": str(INDICES_DIR / "mtb_h37rv"),
    },
    "bcg": {
        "label": "M. bovis BCG Pasteur 1173P2",
        "genome": GENOME_DIR / "bcg_pasteur.fna",
        "index": str(INDICES_DIR / "bcg_pasteur"),
    },
    "msm": {
        "label": "M. smegmatis mc2 155",
        "genome": GENOME_DIR / "msm_mc2155.fna",
        "index": str(INDICES_DIR / "msm_mc2155"),
    },
}
# 背景基因组方向性：设计某物种的探针时，用另外两个物种过滤交叉反应
BACKGROUNDS = {
    "mtb": ["bcg", "msm"],
    "bcg": ["mtb", "msm"],
    "msm": ["mtb", "bcg"],
}


def parse_attributes(field: str) -> dict[str, str]:
    """解析 GFF 第 9 列的 key=value;key=value 属性串。"""
    attrs: dict[str, str] = {}
    for piece in field.strip().rstrip(";").split(";"):
        if "=" in piece:
            key, _, value = piece.partition("=")
            attrs[key.strip()] = value.strip()
    return attrs


def _is_propionyl_gene(gene: str, text: str) -> bool:
    """按 docstring 中的三条标准判断是否为丙酰辅酶A代谢相关基因。"""
    if KEYWORDS.search(text):
        return True
    gene_clean = gene.lower().strip()
    return gene_clean in CORE_GENES


def select_genes_from_mycobrowser(gff_path: Path) -> list[dict]:
    """解析 MycoBrowser GFF（MTB/MSM），返回按 locus 排序的候选基因列表。

    MycoBrowser 的 CDS 属性串形如：
        Locus=Rv1131;Name=prpC;Function=...;Product=...
    """
    genes: list[dict] = []
    with gff_path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attrs = parse_attributes(parts[8])
            locus = attrs.get("Locus", "")
            gene = attrs.get("Name", "") or locus
            # 跳过假基因与 RNA 特征（只设计蛋白编码基因）
            if attrs.get("Is_Pseudogene", "No") == "Yes":
                continue
            text = " ".join(
                attrs.get(key, "") for key in ("Function", "Product", "Comments")
            )
            if not _is_propionyl_gene(gene, text):
                continue
            genes.append(
                {
                    "locus": locus,
                    "gene": gene,
                    "product": attrs.get("Product", "") or attrs.get("Function", ""),
                    "start": int(parts[3]),
                    "stop": int(parts[4]),
                    "strand": parts[6],
                }
            )
    return sorted(genes, key=lambda g: g["locus"])


def select_genes_from_ncbi_gff(gff_path: Path) -> list[dict]:
    """解析 NCBI RefSeq GFF（BCG），筛选标准与 MycoBrowser 相同。

    NCBI 的 CDS 属性串形如：
        ID=cds-BCG_RS06165;Name=prpC;locus_tag=BCG_RS06165;product=...
    """
    genes: list[dict] = []
    seen: set[str] = set()
    with gff_path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attrs = parse_attributes(parts[8])
            locus = attrs.get("locus_tag", "")
            gene = attrs.get("gene", "") or locus
            product = attrs.get("product", "")
            if locus in seen:  # 多片段 CDS（核糖体移码）只取一次
                continue
            seen.add(locus)
            if not _is_propionyl_gene(gene, product):
                continue
            genes.append(
                {
                    "locus": locus,
                    "gene": gene,
                    "product": product,
                    "start": int(parts[3]),
                    "stop": int(parts[4]),
                    "strand": parts[6],
                }
            )
    return sorted(genes, key=lambda g: g["locus"])


def load_genome(path: Path) -> str:
    """读入单记录基因组 FASTA 的序列部分（三个基因组均只有一条染色体）。"""
    pieces: list[str] = []
    with path.open() as fh:
        for line in fh:
            if not line.startswith(">"):
                pieces.append(line.strip())
    return "".join(pieces).upper()


def extract_gene(genome_seq: str, gene: dict) -> str:
    """按 GFF 坐标（1-based 闭区间）提取基因序列；负链取反向互补。"""
    seq = genome_seq[gene["start"] - 1 : gene["stop"]]
    return reverse_complement(seq) if gene["strand"] == "-" else seq


def write_fasta(path: Path, header: str, seq: str) -> None:
    """按每行 80 字符写出 FASTA。"""
    with path.open("w") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i : i + 80] + "\n")


def failure_category(probe) -> str:
    """把一条候选的失败原因归入漏斗阶段（热力学 / 特异性 / 未入选）。"""
    thermo_keys = ("GC=", "Tm=", "homopolymer", "hairpin")
    if any(k in r for r in probe.failure_reasons for k in thermo_keys):
        return "thermo"
    if any("hits" in r for r in probe.failure_reasons):
        return "specificity"
    if any("not_selected" in r for r in probe.failure_reasons):
        return "not_selected"
    return "other"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    GENE_FASTA_DIR.mkdir(parents=True, exist_ok=True)

    # ---- 第 1 步：解析注释、选取基因 ----
    selected: dict[str, list[dict]] = {
        "mtb": select_genes_from_mycobrowser(
            ANNOT_DIR / "mycobrowser_h37rv.gff"
        )[:MAX_PER_ORGANISM],
        "msm": select_genes_from_mycobrowser(
            ANNOT_DIR / "mycobrowser_smegmatis.gff"
        )[:MAX_PER_ORGANISM],
        "bcg": select_genes_from_ncbi_gff(
            ANNOT_DIR / "GCF_000009445.1_ASM944v1_genomic.gff"
        )[:MAX_PER_ORGANISM],
    }
    total = sum(len(v) for v in selected.values())
    print(f"[1/3] 基因选取完成：", flush=True)
    for key, genes in selected.items():
        print(f"      {ORGANISMS[key]['label']}: {len(genes)} 个", flush=True)
    print(f"      合计 {total} 个（要求 ≥30）", flush=True)
    if total < 30:
        raise SystemExit("选取的基因不足 30 个，请检查筛选条件。")

    # ---- 第 2 步：提取基因序列 ----
    genome_seqs = {
        key: load_genome(cfg["genome"]) for key, cfg in ORGANISMS.items()
    }
    gene_meta: list[dict] = []
    for org_key, genes in selected.items():
        for gene in genes:
            seq = extract_gene(genome_seqs[org_key], gene)
            fasta_id = f"{org_key}|{gene['locus']}|{gene['gene']}"
            write_fasta(GENE_FASTA_DIR / f"{org_key}_{gene['locus']}.fasta", fasta_id, seq)
            gene_meta.append(
                {
                    "org": org_key,
                    "locus": gene["locus"],
                    "gene": gene["gene"],
                    "product": gene["product"][:90],
                    "length_nt": len(seq),
                    "fasta_id": fasta_id,
                    "fasta_path": str(GENE_FASTA_DIR / f"{org_key}_{gene['locus']}.fasta"),
                }
            )
    print(f"[2/3] 序列提取完成（{len(gene_meta)} 条），开始批量设计…", flush=True)

    # ---- 第 3 步：批量运行 smFISH 设计 ----
    params = DesignParams(design_scheme="smFISH", desired_probe_count=20)
    rows: list[dict] = []
    results_by_gene: dict[str, tuple[DesignResult, list[ReferenceGenome]]] = {}
    started = time.time()
    for i, meta in enumerate(gene_meta, 1):
        org_key = meta["org"]
        target_cfg = ORGANISMS[org_key]
        # 背景基因组 = 另外两个物种
        hosts = [
            ReferenceGenome(
                id=bg_key,
                organism=ORGANISMS[bg_key]["label"],
                fasta_path=str(ORGANISMS[bg_key]["genome"]),
                bowtie2_index=ORGANISMS[bg_key]["index"],
                is_host=True,
            )
            for bg_key in BACKGROUNDS[org_key]
        ]
        t0 = time.time()
        try:
            result: DesignResult = run_design(
                meta["fasta_path"], target_cfg["index"], hosts, params, threads=4
            )
        except Exception as exc:  # 记录失败基因，继续批量
            rows.append({**meta, "error": str(exc)[:200]})
            print(f"  ({i}/{len(gene_meta)}) {org_key}:{meta['gene']} 失败：{exc}", flush=True)
            continue
        results_by_gene[meta["fasta_id"]] = (result, hosts)

        passed = result.passed_probes
        tms = [p.tm for p in passed] or [0.0]
        # 漏斗统计：每阶段的存活候选数
        thermo_pass = sum(1 for p in result.probes if failure_category(p) not in ("thermo",))
        spec_pass = sum(
            1 for p in result.probes
            if failure_category(p) not in ("thermo", "specificity")
        )
        # 背景淘汰专列：统计被宿主/背景基因组淘汰的候选数（工作流测试的重点观察项）
        host_rejected = sum(
            1 for p in result.probes
            if any("host_hits" in r for r in p.failure_reasons)
        )
        target_rejected = sum(
            1 for p in result.probes
            if any("target_hits" in r for r in p.failure_reasons)
        )
        row = {
            **meta,
            "candidates": len(result.probes),
            "thermo_pass": thermo_pass,
            "spec_pass": spec_pass,
            "final_probes": len(passed),
            "coverage_pct": round(
                sum(p.stop - p.start for p in passed) / result.target_length * 100, 1
            ) if result.target_length else 0.0,
            "tm_mean": round(sum(tms) / len(tms), 1) if passed else "",
            "host_rejected": host_rejected,
            "target_rejected": target_rejected,
            "runtime_s": round(time.time() - t0, 1),
            "error": "",
        }
        rows.append(row)
        print(
            f"  ({i}/{len(gene_meta)}) {org_key}:{meta['gene']} "
            f"候选{row['candidates']} → 最终{row['final_probes']} "
            f"(背景淘汰{host_rejected}, {row['runtime_s']}s)",
            flush=True,
        )

    # ---- 第 4 步：落盘汇总 ----
    with (OUT_DIR / "batch_summary.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    # 每个基因的最终探针明细（含序列与坐标），合并为一张总表（复用第 3 步结果）
    detail_rows: list[dict] = []
    for meta, row in zip(gene_meta, rows):
        if row.get("error") or meta["fasta_id"] not in results_by_gene:
            continue
        result, _hosts = results_by_gene[meta["fasta_id"]]
        for probe in result.passed_probes:
            detail_rows.append(
                {
                    "org": meta["org"],
                    "locus": meta["locus"],
                    "gene": meta["gene"],
                    "probe_id": probe.probe_id,
                    "start": probe.start,
                    "stop": probe.stop,
                    "sequence": probe.sequence,
                    "tm": round(probe.tm, 1),
                    "gc": round(probe.gc_content, 2),
                    "target_hits": probe.target_hits,
                    "host_hits": json.dumps(probe.host_hits),
                    "score": round(probe.score, 3),
                }
            )
    with (OUT_DIR / "all_final_probes.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(detail_rows[0].keys()))
        writer.writeheader()
        writer.writerows(detail_rows)

    elapsed = time.time() - started
    print(f"[3/3] 完成：{len(gene_meta)} 个基因，总耗时 {elapsed / 60:.1f} 分钟", flush=True)
    print(f"      汇总：{OUT_DIR / 'batch_summary.csv'}", flush=True)
    print(f"      探针明细：{OUT_DIR / 'all_final_probes.csv'}", flush=True)


if __name__ == "__main__":
    main()
