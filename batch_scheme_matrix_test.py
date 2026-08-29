#!/usr/bin/env python
"""在 48 个低丰度基因上测试 smiFISH / HCR 3.0 / SNAIL 方案的可行性。

复用第二轮测试（batch_lowabundance_test.py）选取的低丰度基因集与其输入
FASTA（test_data/lowabundance_results/gene_fastas/），对每个基因分别以
smiFISH、HCR 3.0、SNAIL FISH 三种方案各跑一次完整设计，回答"每种方案
在这些基因上是否可行、产出多少、有什么限制"。

各方案配置（与实验目的匹配）：
    smiFISH   与 smFISH 相同的堆叠逻辑（低丰度 → 目标 48 条），readout
              使用占位序列（订购前替换为实验所用 LNA 二级探针的互补序列）；
    HCR 3.0   信号由放大器链扩增，每转录本 20 对半探针即可，initiator 通道
              B1，GC 45–55%、Gibbs −70~−50 kcal/mol、dTm ≤ 5 °C（协议默认）；
    SNAIL     连接+滚环扩增，20 对 primer/padlock；UGI 条码留空（N 占位），
              padlock 订购需加 /5Phos/。

同时做方案产物的结构自检：
    smiFISH  full_sequence == 探针 + linker + readout（3′ 方案）
    HCR3     P1 = initiator(odd) + 3′半探针，P2 = 5′半探针 + initiator(even)
    SNAIL    primer 以臂1开头；padlock 以 5′ anchor 开头、含 UGI 占位

运行（项目目录）：
    /opt/anaconda3/envs/Probe/bin/python batch_scheme_matrix_test.py
"""

from __future__ import annotations

import csv
import statistics
import sys
import time
from pathlib import Path

PROJECT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT / "src"))

from batch_lowabundance_test import ORGANISMS  # noqa: E402
from batch_smfish_test import failure_category  # noqa: E402
from probedesign.models import DesignParams, DesignResult  # noqa: E402
from probedesign.pipeline import run_design  # noqa: E402
from probedesign.schemes.initiators import HCR_INITIATORS  # noqa: E402

GENE_CSV = PROJECT / "test_data" / "lowabundance_results" / "batch_summary.csv"
OUT_DIR = PROJECT / "test_data" / "scheme_matrix_results"

# smiFISH 占位 readout：仅用于工作流测试，订购前替换为实际 LNA 二级探针互补序列
PLACEHOLDER_READOUT = "ACGTCGACTATCGAT"

SCHEMES: dict[str, dict] = {
    "smiFISH": {
        "label": "smiFISH",
        "rationale": "低丰度堆叠逻辑与 smFISH 相同 → 目标 48 条/基因",
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
        "label": "HCR 3.0",
        "rationale": "放大器扩增信号 → 20 对半探针/基因即可",
        "params": lambda: DesignParams(
            design_scheme="HCR3",
            desired_probe_count=20,
            min_gap=2,
            hcr_channel="B1",
        ),
    },
    "SNAIL-FISH": {
        "label": "SNAIL FISH",
        "rationale": "连接+滚环扩增 → 20 对 primer/padlock/基因；UGI 用 N 占位",
        "params": lambda: DesignParams(
            design_scheme="SNAIL-FISH",
            desired_probe_count=20,
            min_gap=2,
        ),
    },
}


def scheme_sanity_check(scheme: str, result: DesignResult) -> list[str]:
    """对通过探针做产物结构自检，返回问题列表（空列表 = 全部通过）。"""
    problems: list[str] = []
    params = result.params
    for probe in result.passed_probes:
        if scheme == "smiFISH":
            full = probe.metadata.get("full_sequence", "")
            expected_tail = probe.sequence + "TTT" + PLACEHOLDER_READOUT
            if full != expected_tail:
                problems.append(f"{probe.probe_id}: full_sequence 与预期拼接不一致")
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
                problems.append(f"{probe.probe_id}: dTm={dtm} 超过上限却在最终集合")
        elif scheme == "SNAIL-FISH":
            primer = probe.metadata.get("primer_sequence", "")
            padlock = probe.metadata.get("padlock_sequence", "")
            if not primer.startswith(probe.metadata.get("arm1_sequence", "??")):
                problems.append(f"{probe.probe_id}: primer 未以臂1开头")
            if not padlock.startswith(params.snail_padlock_start):
                problems.append(f"{probe.probe_id}: padlock 未以 5′ anchor 开头")
            if "N" not in padlock:
                problems.append(f"{probe.probe_id}: padlock 缺少 UGI 占位")
    return problems


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    genes = list(csv.DictReader(open(GENE_CSV)))
    genes = [g for g in genes if not g.get("error")]
    print(f"载入 {len(genes)} 个低丰度基因，测试 {len(SCHEMES)} 个方案…", flush=True)

    rows: list[dict] = []
    started = time.time()
    for scheme_key, cfg in SCHEMES.items():
        params = cfg["params"]()
        for i, gene in enumerate(genes, 1):
            t0 = time.time()
            try:
                result: DesignResult = run_design(
                    gene["fasta_path"],
                    ORGANISMS[gene["org"]]["index"],
                    [],  # 低丰度独立检测：无跨物种背景
                    params,
                    threads=4,
                )
            except Exception as exc:
                rows.append(
                    {"scheme": scheme_key, **gene, "error": str(exc)[:200]}
                )
                print(f"  [{scheme_key}] ({i}) {gene['org']}:{gene['gene']} 失败：{exc}", flush=True)
                continue
            passed = result.passed_probes
            tms = [p.tm for p in passed]
            problems = scheme_sanity_check(scheme_key, result)
            rows.append(
                {
                    "scheme": scheme_key,
                    "org": gene["org"],
                    "gene": gene["gene"],
                    "locus": gene["locus"],
                    "cls": gene["cls"],
                    "length_nt": gene["length_nt"],
                    "candidates": len(result.probes),
                    "thermo_pass": sum(
                        1 for p in result.probes
                        if failure_category(p) != "thermo"
                    ),
                    "final_probes": len(passed),
                    "coverage_pct": round(
                        sum(p.stop - p.start for p in passed) / result.target_length * 100, 1
                    ) if result.target_length else 0.0,
                    "tm_mean": round(statistics.mean(tms), 1) if tms else "",
                    "sanity_problems": "; ".join(problems)[:200],
                    "runtime_s": round(time.time() - t0, 1),
                    "error": "",
                }
            )
        print(f"[{cfg['label']}] 48 基因完成（累计 {time.time() - started:.0f}s）", flush=True)

    with (OUT_DIR / "scheme_matrix_summary.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    elapsed = time.time() - started
    print(f"完成：{len(SCHEMES)} 方案 × {len(genes)} 基因，总耗时 {elapsed:.0f} 秒", flush=True)
    print(f"汇总：{OUT_DIR / 'scheme_matrix_summary.csv'}", flush=True)


if __name__ == "__main__":
    main()
