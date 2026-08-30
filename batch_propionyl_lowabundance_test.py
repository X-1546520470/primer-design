#!/usr/bin/env python
"""丙酰辅酶A代谢基因 · 低丰度转录本检测场景下的四方案批量测试。

实验目的（与第一/二轮物种区分测试不同）：
    **液培、低菌量条件下，独立检测各分枝杆菌（MTB / BCG / MSM）低丰度靶
    基因的转录表达。** 每种菌独立培养、独立制片检测，样本中只存在该菌种
    自身的 RNA，不存在另一个物种的核酸，因此：
        - **关闭跨物种背景过滤**（host_genomes = []）——不需要区分物种，
          也不应该用另一个物种的基因组淘汰探针（否则保守基因会被误杀）；
        - 保留**靶标自身基因组的重复序列 QC**（max_target_hits=10），这是
          探针质量的内在检查，与跨物种过滤无关。

低丰度检测的参数取向（信号 ∝ 结合在单个转录本上的探针数）：
    - smFISH / smiFISH：目标 48 条/转录本（经典 smFISH 低丰度堆叠），
      min_gap=2；
    - HCR 3.0：信号由放大器链指数扩增，20 对/转录本即可，min_gap=2；
    - SNAIL：连接+滚环扩增，20 对 primer/padlock/转录本，min_gap=2。

基因集：与第一轮相同的丙酰辅酶A代谢相关基因（MycoBrowser + NCBI 注释，
筛选标准见 batch_smfish_test.py），重新选取并提取序列，共 40 个
（MTB 16 + MSM 15 + BCG 9，≥30）。

产物：每个 (方案 × 基因) 记录候选/热力学存活/特异性存活/最终条数/覆盖率/
Tm 均一性/靶重复淘汰/运行时长/异常，并做方案专属结构自检与警告捕获。
本脚本只测试、不修改引擎。

运行（项目目录）：
    /opt/anaconda3/envs/Probe/bin/python batch_propionyl_lowabundance_test.py
"""

from __future__ import annotations

import csv
import statistics
import sys
import time
import warnings
from pathlib import Path

PROJECT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT / "src"))

# 复用第一轮的基因选取 / 序列提取 / 物种配置工具
from batch_smfish_test import (  # noqa: E402
    ANNOT_DIR,
    ORGANISMS,
    extract_gene,
    failure_category,
    load_genome,
    select_genes_from_mycobrowser,
    select_genes_from_ncbi_gff,
    write_fasta,
)
from mycoprimer.models import DesignParams, DesignResult  # noqa: E402
from mycoprimer.pipeline import run_design  # noqa: E402
from mycoprimer.schemes.initiators import HCR_INITIATORS  # noqa: E402

OUT_DIR = PROJECT / "test_data" / "propionyl_lowabundance_results"
GENE_FASTA_DIR = OUT_DIR / "gene_fastas"
MAX_PER_ORGANISM = 16

# smiFISH 占位 readout（工作流测试用；订购前替换为实际 LNA 二级探针互补序列）
PLACEHOLDER_READOUT = "ACGTCGACTATCGAT"

# 低丰度检测配置：无跨物种背景；线性探针高堆叠(48)，扩增方案 20 对
SCHEMES: dict[str, dict] = {
    "smFISH": {
        "label": "smFISH",
        "params": lambda: DesignParams(
            design_scheme="smFISH",
            desired_probe_count=48,
            min_gap=2,
        ),
    },
    "smiFISH": {
        "label": "smiFISH",
        "params": lambda: DesignParams(
            design_scheme="smiFISH",
            desired_probe_count=48,
            min_gap=2,
            smi_readout_sequence=PLACEHOLDER_READOUT,
            smi_readout_position="3prime",
            smi_linker="TTT",
        ),
    },
    "HCR3": {
        "label": "HCR3",
        "params": lambda: DesignParams(
            design_scheme="HCR3",
            desired_probe_count=20,
            min_gap=2,
            hcr_channel="B1",
        ),
    },
    "SNAIL-FISH": {
        "label": "SNAIL-FISH",
        "params": lambda: DesignParams(
            design_scheme="SNAIL-FISH",
            desired_probe_count=20,
            min_gap=2,
        ),
    },
}


def scheme_sanity_check(scheme: str, result: DesignResult) -> list[str]:
    """对通过探针做方案专属结构自检，返回问题描述列表（空 = 全部通过）。"""
    problems: list[str] = []
    params = result.params
    n_checked = 0
    for probe in result.passed_probes:
        n_checked += 1
        if scheme == "smiFISH":
            full = probe.metadata.get("full_sequence", "")
            expected = probe.sequence + "TTT" + PLACEHOLDER_READOUT
            if full != expected:
                problems.append(f"{probe.probe_id}: full_sequence 拼接不一致")
            if any(c not in "ACGT" for c in full):
                problems.append(f"{probe.probe_id}: full_sequence 含非 ACGT 字符")
        elif scheme == "HCR3":
            p1 = probe.metadata.get("P1_sequence", "")
            p2 = probe.metadata.get("P2_sequence", "")
            initiator = HCR_INITIATORS[params.hcr_channel]
            if not p1.startswith(initiator["odd"]):
                problems.append(f"{probe.probe_id}: P1 未以 odd initiator 开头")
            if not p2.endswith(initiator["even"]):
                problems.append(f"{probe.probe_id}: P2 未以 even initiator 结尾")
            dtm = probe.metadata.get("dTm")
            if dtm is not None and params.hcr_dtm_max is not None and dtm > params.hcr_dtm_max:
                problems.append(f"{probe.probe_id}: dTm={dtm:.1f} 超上限却进入最终集合")
            if len(p1) - len(initiator["odd"]) != 25 or len(p2) - len(initiator["even"]) != 25:
                problems.append(f"{probe.probe_id}: 半探针靶结合长度不是 25 nt")
        elif scheme == "SNAIL-FISH":
            primer = probe.metadata.get("primer_sequence", "")
            padlock = probe.metadata.get("padlock_sequence", "")
            arm1 = probe.metadata.get("arm1_sequence", "")
            arm2 = probe.metadata.get("arm2_sequence", "")
            if not primer.startswith(arm1 or "??"):
                problems.append(f"{probe.probe_id}: primer 未以臂1开头")
            if not primer.endswith(params.snail_primer_end):
                problems.append(f"{probe.probe_id}: primer 未以 3′ linker 结尾")
            if not padlock.startswith(params.snail_padlock_start):
                problems.append(f"{probe.probe_id}: padlock 未以 5′ anchor 开头")
            if not padlock.endswith(params.snail_padlock_end):
                problems.append(f"{probe.probe_id}: padlock 未以 3′ anchor 结尾")
            if "N" not in padlock:
                problems.append(f"{probe.probe_id}: padlock 缺少 UGI 占位")
            if len(arm1) != params.snail_arm_length or len(arm2) != params.snail_arm_length:
                problems.append(f"{probe.probe_id}: 臂长不是 {params.snail_arm_length} nt")
    if n_checked == 0:
        problems.append("零产出：无通过探针可检")
    return problems


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    GENE_FASTA_DIR.mkdir(parents=True, exist_ok=True)

    # ---- 第 1 步：重新选取基因 ----
    selected = {
        "mtb": select_genes_from_mycobrowser(ANNOT_DIR / "mycobrowser_h37rv.gff")[:MAX_PER_ORGANISM],
        "msm": select_genes_from_mycobrowser(ANNOT_DIR / "mycobrowser_smegmatis.gff")[:MAX_PER_ORGANISM],
        "bcg": select_genes_from_ncbi_gff(ANNOT_DIR / "GCF_000009445.1_ASM944v1_genomic.gff")[:MAX_PER_ORGANISM],
    }
    total = sum(len(v) for v in selected.values())
    print(f"[1/3] 基因选取完成（丙酰辅酶A代谢，低丰度检测场景）：", flush=True)
    for key, genes in selected.items():
        print(f"      {ORGANISMS[key]['label']}: {len(genes)} 个", flush=True)
    print(f"      合计 {total} 个（要求 ≥30）", flush=True)
    if total < 30:
        raise SystemExit("选取的基因不足 30 个，请检查筛选条件。")

    # ---- 第 2 步：提取序列 ----
    genome_seqs = {key: load_genome(cfg["genome"]) for key, cfg in ORGANISMS.items()}
    gene_meta: list[dict] = []
    for org_key, genes in selected.items():
        for gene in genes:
            seq = extract_gene(genome_seqs[org_key], gene)
            fasta_id = f"{org_key}|{gene['locus']}|{gene['gene']}"
            path = GENE_FASTA_DIR / f"{org_key}_{gene['locus']}.fasta"
            write_fasta(path, fasta_id, seq)
            gene_meta.append(
                {
                    "org": org_key,
                    "locus": gene["locus"],
                    "gene": gene["gene"],
                    "product": gene["product"][:90],
                    "length_nt": len(seq),
                    "fasta_path": str(path),
                }
            )
    print(f"[2/3] 序列提取完成（{len(gene_meta)} 条），开始四方案批量设计…", flush=True)

    # ---- 第 3 步：四方案批量设计（无跨物种背景；只做自身基因组重复检测）----
    rows: list[dict] = []
    started = time.time()
    for scheme_key, cfg in SCHEMES.items():
        scheme_started = time.time()
        for i, meta in enumerate(gene_meta, 1):
            org_key = meta["org"]
            params = cfg["params"]()
            t0 = time.time()
            warn_msgs: list[str] = []
            try:
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    result: DesignResult = run_design(
                        meta["fasta_path"],
                        ORGANISMS[org_key]["index"],
                        [],  # 低丰度独立检测：无跨物种背景
                        params,
                        threads=4,
                    )
                for w in caught:
                    warn_msgs.append(f"{w.category.__name__}: {str(w.message)[:80]}")
            except Exception as exc:
                rows.append(
                    {
                        "scheme": scheme_key, "org": org_key, "locus": meta["locus"],
                        "gene": meta["gene"], "length_nt": meta["length_nt"],
                        "candidates": "", "thermo_pass": "", "spec_pass": "",
                        "final_probes": 0, "coverage_pct": "", "tm_mean": "", "tm_sd": "",
                        "target_rejected": "", "zero_output": "YES",
                        "sanity_problems": "run_failed", "warnings": "",
                        "runtime_s": round(time.time() - t0, 1), "error": str(exc)[:300],
                    }
                )
                print(f"  [{scheme_key}] ({i}/{len(gene_meta)}) {org_key}:{meta['gene']} 失败：{exc}", flush=True)
                continue

            thermo_pass = sum(1 for p in result.probes if failure_category(p) != "thermo")
            spec_pass = sum(1 for p in result.probes if failure_category(p) not in ("thermo", "specificity"))
            target_rejected = sum(
                1 for p in result.probes if any("target_hits" in r for r in p.failure_reasons)
            )
            passed = result.passed_probes
            tms = [p.tm for p in passed]
            problems = scheme_sanity_check(scheme_key, result)
            coverage = (
                round(sum(p.stop - p.start for p in passed) / result.target_length * 100, 1)
                if result.target_length else 0.0
            )
            rows.append(
                {
                    "scheme": scheme_key, "org": org_key, "locus": meta["locus"],
                    "gene": meta["gene"], "length_nt": meta["length_nt"],
                    "candidates": len(result.probes),
                    "thermo_pass": thermo_pass, "spec_pass": spec_pass,
                    "final_probes": len(passed), "coverage_pct": coverage,
                    "tm_mean": round(statistics.mean(tms), 1) if tms else "",
                    "tm_sd": round(statistics.stdev(tms), 1) if len(tms) > 1 else 0.0,
                    "target_rejected": target_rejected,
                    "zero_output": "YES" if len(passed) == 0 else "",
                    "sanity_problems": "; ".join(problems)[:300],
                    "warnings": " | ".join(sorted(set(warn_msgs)))[:200],
                    "runtime_s": round(time.time() - t0, 1), "error": "",
                }
            )
        print(
            f"[{cfg['label']}] {len(gene_meta)} 基因完成"
            f"（方案耗时 {time.time() - scheme_started:.0f}s，累计 {time.time() - started:.0f}s）",
            flush=True,
        )

    out_csv = OUT_DIR / "propionyl_lowabundance_summary.csv"
    with out_csv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    elapsed = time.time() - started
    n_fail = sum(1 for r in rows if r["error"])
    n_zero = sum(1 for r in rows if r["zero_output"] == "YES" and not r["error"])
    n_sanity = sum(1 for r in rows if r["sanity_problems"] and r["sanity_problems"] != "run_failed")
    print("=" * 60, flush=True)
    print(f"完成：{len(rows)} 个 (方案×基因) 运行，总耗时 {elapsed / 60:.1f} 分钟", flush=True)
    print(f"  运行异常：{n_fail}；零产出：{n_zero}；结构自检报问题：{n_sanity}", flush=True)
    print(f"汇总：{out_csv}", flush=True)


if __name__ == "__main__":
    main()
