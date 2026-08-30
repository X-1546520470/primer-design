#!/usr/bin/env python
"""丙酰辅酶A（propionyl-CoA）代谢基因集上的四方案探针设计批量测试。

测试对象：第一轮（batch_smfish_test.py）从 MycoBrowser / NCBI 注释中选取的
40 个丙酰辅酶A代谢相关基因（MTB H37Rv 16 + M. smegmatis 15 + BCG Pasteur 9），
序列 FASTA 已在 test_data/batch_results/gene_fastas/ 下。本脚本复用该基因集，
对四套工作流各跑一次完整设计：

    smFISH    18–24 nt 反义寡核苷酸，目标 20 条/基因
    smiFISH   smFISH + 共享 readout 延伸段（占位 readout），目标 20 条/基因
    HCR3      52 nt tile 拆两条 25-mer 半探针 + 分裂 initiator，目标 20 对/基因
    SNAIL-FISH primer + 5′ 磷酸化 padlock 双臂结构，目标 20 对/基因

背景基因组设置（与第一轮 smFISH 测试一致，考察交叉反应过滤）：
    MTB 基因 → 背景 = BCG + MSM
    BCG 基因 → 背景 = MTB + MSM
    MSM 基因 → 背景 = MTB + BCG

每个 (方案 × 基因) 记录：候选数、热力学存活、特异性存活、最终条数、覆盖率、
Tm 均一性、背景/靶标淘汰数、运行时长、异常信息；并对通过探针做方案专属
结构自检（smiFISH 拼接、HCR3 initiator/dTm、SNAIL primer/padlock 结构）。
运行期警告（warnings）一并捕获。

注意：本脚本只做测试与记录，不修改引擎；发现的系统性问题在报告中列出。

运行（项目目录）：
    /opt/anaconda3/envs/Probe/bin/python batch_propionyl_matrix_test.py
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

from batch_smfish_test import ORGANISMS, BACKGROUNDS  # noqa: E402
from mycoprimer.models import DesignParams, DesignResult, ReferenceGenome  # noqa: E402
from mycoprimer.pipeline import run_design  # noqa: E402
from mycoprimer.schemes.initiators import HCR_INITIATORS  # noqa: E402

GENE_CSV = PROJECT / "test_data" / "batch_results" / "batch_summary.csv"
OUT_DIR = PROJECT / "test_data" / "propionyl_matrix_results"

# smiFISH 占位 readout（工作流测试用；订购前替换为实际 LNA 二级探针互补序列）
PLACEHOLDER_READOUT = "ACGTCGACTATCGAT"

SCHEMES: dict[str, dict] = {
    "smFISH": {
        "label": "smFISH",
        "params": lambda: DesignParams(
            design_scheme="smFISH",
            desired_probe_count=20,
            min_gap=2,
        ),
    },
    "smiFISH": {
        "label": "smiFISH",
        "params": lambda: DesignParams(
            design_scheme="smiFISH",
            desired_probe_count=20,
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


def failure_stage(probe) -> str:
    """把一条候选归入最终失败阶段：thermo / specificity / not_selected / other / pass。

    用 failure_reasons 关键词分类，覆盖四套方案的原因串：
        热力学  GC= / Tm= / homopolymer / hairpin / Gibbs / dTm / arm1|arm2
        特异性  hits（target_hits / host_hits / primer_host_hits 等）
        选点    not_selected
    """
    if probe.passed:
        return "pass"
    reasons = " ".join(probe.failure_reasons)
    thermo_keys = ("GC=", "Tm=", "homopolymer", "hairpin", "Gibbs", "dTm", "arm1", "arm2")
    if any(k in reasons for k in thermo_keys):
        return "thermo"
    if "hits" in reasons:
        return "specificity"
    if "not_selected" in reasons:
        return "not_selected"
    return "other"


def scheme_sanity_check(scheme: str, result: DesignResult) -> list[str]:
    """对通过探针做方案专属结构自检，返回问题描述列表（空 = 全部通过）。"""
    problems: list[str] = []
    params = result.params
    n_checked = 0
    for probe in result.passed_probes:
        n_checked += 1
        if scheme == "smiFISH":
            full = probe.metadata.get("full_sequence", "")
            expected_tail = probe.sequence + "TTT" + PLACEHOLDER_READOUT
            if full != expected_tail:
                problems.append(f"{probe.probe_id}: full_sequence 与 探针+linker+readout 拼接不一致")
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
            # 半探针长度应为 25 nt（initiator 之外的靶结合部分）
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

    genes = list(csv.DictReader(GENE_CSV.open()))
    genes = [g for g in genes if not g.get("error")]
    # 校验 FASTA 存在
    missing = [g["fasta_id"] for g in genes if not Path(g["fasta_path"]).is_file()]
    if missing:
        raise SystemExit(f"缺少基因 FASTA：{missing[:5]} …")
    print(f"载入 {len(genes)} 个丙酰辅酶A基因，测试 {len(SCHEMES)} 套工作流…", flush=True)

    rows: list[dict] = []
    started = time.time()
    for scheme_key, cfg in SCHEMES.items():
        scheme_started = time.time()
        for i, gene in enumerate(genes, 1):
            org_key = gene["org"]
            hosts = [
                ReferenceGenome(
                    id=bg,
                    organism=ORGANISMS[bg]["label"],
                    fasta_path=str(ORGANISMS[bg]["genome"]),
                    bowtie2_index=ORGANISMS[bg]["index"],
                    is_host=True,
                )
                for bg in BACKGROUNDS[org_key]
            ]
            params = cfg["params"]()
            t0 = time.time()
            warn_msgs: list[str] = []
            try:
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    result: DesignResult = run_design(
                        gene["fasta_path"],
                        ORGANISMS[org_key]["index"],
                        hosts,
                        params,
                        threads=4,
                    )
                for w in caught:
                    warn_msgs.append(f"{w.category.__name__}: {str(w.message)[:80]}")
            except Exception as exc:  # 记录失败，继续批量
                rows.append(
                    {
                        "scheme": scheme_key,
                        "org": org_key,
                        "locus": gene["locus"],
                        "gene": gene["gene"],
                        "length_nt": gene["length_nt"],
                        "candidates": "",
                        "thermo_pass": "",
                        "spec_pass": "",
                        "final_probes": 0,
                        "coverage_pct": "",
                        "tm_mean": "",
                        "tm_sd": "",
                        "host_rejected": "",
                        "target_rejected": "",
                        "zero_output": "YES",
                        "sanity_problems": "run_failed",
                        "warnings": "",
                        "runtime_s": round(time.time() - t0, 1),
                        "error": str(exc)[:300],
                    }
                )
                print(
                    f"  [{scheme_key}] ({i}/{len(genes)}) {org_key}:{gene['gene']} 运行失败：{exc}",
                    flush=True,
                )
                continue

            stages = [failure_stage(p) for p in result.probes]
            thermo_pass = sum(1 for s in stages if s not in ("thermo",))
            spec_pass = sum(1 for s in stages if s not in ("thermo", "specificity"))
            host_rejected = sum(
                1 for p in result.probes
                if any("host_hits" in r for r in p.failure_reasons)
            )
            target_rejected = sum(
                1 for p in result.probes
                if any("target_hits" in r for r in p.failure_reasons)
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
                    "scheme": scheme_key,
                    "org": org_key,
                    "locus": gene["locus"],
                    "gene": gene["gene"],
                    "length_nt": gene["length_nt"],
                    "candidates": len(result.probes),
                    "thermo_pass": thermo_pass,
                    "spec_pass": spec_pass,
                    "final_probes": len(passed),
                    "coverage_pct": coverage,
                    "tm_mean": round(statistics.mean(tms), 1) if tms else "",
                    "tm_sd": round(statistics.stdev(tms), 1) if len(tms) > 1 else 0.0,
                    "host_rejected": host_rejected,
                    "target_rejected": target_rejected,
                    "zero_output": "YES" if len(passed) == 0 else "",
                    "sanity_problems": "; ".join(problems)[:300],
                    "warnings": " | ".join(sorted(set(warn_msgs)))[:200],
                    "runtime_s": round(time.time() - t0, 1),
                    "error": "",
                }
            )
        print(
            f"[{cfg['label']}] {len(genes)} 基因完成"
            f"（方案耗时 {time.time() - scheme_started:.0f}s，累计 {time.time() - started:.0f}s）",
            flush=True,
        )

    out_csv = OUT_DIR / "propionyl_matrix_summary.csv"
    with out_csv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    elapsed = time.time() - started
    n_runs = len(rows)
    n_fail = sum(1 for r in rows if r["error"])
    n_zero = sum(1 for r in rows if r["zero_output"] == "YES" and not r["error"])
    n_sanity = sum(1 for r in rows if r["sanity_problems"] and r["sanity_problems"] != "run_failed")
    print("=" * 60, flush=True)
    print(f"完成：{n_runs} 个 (方案×基因) 运行，总耗时 {elapsed / 60:.1f} 分钟", flush=True)
    print(f"  运行异常：{n_fail}；零产出（无异常但无探针）：{n_zero}；结构自检报问题：{n_sanity}", flush=True)
    print(f"汇总：{out_csv}", flush=True)


if __name__ == "__main__":
    main()
