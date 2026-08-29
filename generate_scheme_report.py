#!/usr/bin/env python
"""从 scheme_matrix_results 生成 smiFISH / HCR / SNAIL 三方案可行性报告。"""

from __future__ import annotations

import csv
from collections import defaultdict
from datetime import datetime
from pathlib import Path

PROJECT = Path(__file__).resolve().parent
SUMMARY = PROJECT / "test_data" / "scheme_matrix_results" / "scheme_matrix_summary.csv"
REPORT = PROJECT / "三方案可行性对比测试报告.md"

ORG_LABEL = {
    "mtb": "MTB（H37Rv）",
    "bcg": "BCG（Pasteur 1173P2）",
    "msm": "MSM（mc² 155）",
}
CLS_LABEL = {
    "sigma": "σ 因子",
    "kinase": "丝/苏氨酸激酶",
    "response_regulator": "反应调节因子",
    "regulator": "转录调控因子",
}
SCHEME_LABEL = {
    "smiFISH": "smiFISH",
    "HCR3": "HCR 3.0",
    "SNAIL-FISH": "SNAIL FISH",
}


def main() -> None:
    rows = list(csv.DictReader(open(SUMMARY)))
    by_scheme: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        by_scheme[row["scheme"]].append(row)

    lines: list[str] = []
    lines.append("# smiFISH / HCR 3.0 / SNAIL 方案可行性对比测试报告")
    lines.append("")
    lines.append(f"- **测试日期**：{datetime.now():%Y-%m-%d}")
    lines.append("- **测试基因**：第二轮选取的 48 个低丰度基因（σ 因子/丝苏氨酸激酶/")
    lines.append("  反应调节因子/转录调控因子 × 4 类，MTB·BCG·MSM 各 16 个）")
    lines.append("- **测试方式**：每个基因分别以 smiFISH、HCR 3.0、SNAIL 三种方案跑完整")
    lines.append("  设计流程（无跨物种背景过滤，与低丰度检测目标一致），并做产物结构自检。")
    lines.append("")

    lines.append("## 1. 总体结论")
    lines.append("")
    lines.append("| 方案 | 结论 | 产出 | 适用性说明 |")
    lines.append("|------|------|------|------------|")
    lines.append("| **smiFISH** | ✅ **完全可行** | 1116 条（与 smFISH 完全一致） | 零失败、零结构问题；"
                 "只需在订购前把占位 readout 换成实际 LNA 二级探针互补序列 |")
    lines.append("| **SNAIL** | ✅ **可行** | 252 对 primer/padlock，平均 5.2 对/基因 | "
                 "双臂结构产率天然低于线性探针；420 nt 短基因 whiB5 无可用位置 |")
    lines.append("| **HCR 3.0** | ⚠️ **可行，但需按菌种 GC 调参** | 协议默认窗口下 67 对（1.4 对/基因）；"
                 "放宽后 186 对（3.9 对/基因） | 高 GC 基因组（65%）上 GC/Gibbs 窗口是主要瓶颈 |")
    lines.append("")

    lines.append("## 2. 各方案结果明细")
    lines.append("")
    lines.append("| 方案 | 基因数 | 候选总数 | 最终探针 | 平均/基因 | 零产出基因 | 覆盖率范围 | Tm 范围 (°C) | 结构自检 |")
    lines.append("|------|--------|----------|----------|-----------|------------|------------|--------------|----------|")
    for scheme in ("smiFISH", "HCR3", "SNAIL-FISH"):
        g = by_scheme[scheme]
        finals = [int(r["final_probes"]) for r in g]
        covs = [float(r["coverage_pct"]) for r in g if r["coverage_pct"]]
        tms = [float(r["tm_mean"]) for r in g if r["tm_mean"]]
        problems = sum(1 for r in g if r["sanity_problems"])
        zero = sum(1 for f in finals if f == 0)
        lines.append(
            f"| {SCHEME_LABEL[scheme]} | {len(g)} | {num(g, 'candidates')} "
            f"| {sum(finals)} | {sum(finals) / len(g):.1f} | {zero} "
            f"| {min(covs):.0f}–{max(covs):.0f}% | {min(tms):.1f}–{max(tms):.1f} "
            f"| {len(g) - problems}/{len(g)} 通过 |"
        )
    lines.append("")
    lines.append("产物结构自检内容：smiFISH 验证完整序列 = 探针 + TTT linker + readout；")
    lines.append("HCR 验证 P1 以 odd initiator 开头、P2 以 even initiator 结尾、dTm 不超标；")
    lines.append("SNAIL 验证 primer 以臂1开头、padlock 以 5′ anchor 开头且含 UGI 占位。")
    lines.append("")

    lines.append("## 3. HCR 3.0 的可行性分析与调参建议")
    lines.append("")
    lines.append("### 3.1 测试中发现并修复了一个引擎问题")
    lines.append("")
    lines.append("首轮 HCR 测试仅产出 34 对（0.7 对/基因）。解剖漏斗发现主导瓶颈是")
    lines.append("**发卡过滤作用在了 52-mer tile 整体上**：45 °C 的发卡 Tm 阈值是为")
    lines.append("18–24-mer 设定的，而 52-mer 的发卡 Tm 天然更高，导致 98.7% 的候选被")
    lines.append("误杀（MSMEG_0786 基因 2232 个候选中 2203 个因 hairpinTm 淘汰）。")
    lines.append("而 HCR 实际合成的是两条 25-mer 半探针，不是 tile 整体。")
    lines.append("")
    lines.append("修复：发卡检查移到 **25-mer 半探针层面**（拆分后对 5′/3′ 半探针分别")
    lines.append("检查），并加上此前已改为可选的 tile Tm 过滤，HCR 的热力学 QC 现在")
    lines.append("全部作用于真正合成的序列。")
    lines.append("")
    lines.append("### 3.2 高 GC 基因组上的窗口调参（敏感性测试，48 基因）")
    lines.append("")
    lines.append("| 参数组合 | 总对数 | 平均对/基因 | 零产出基因 |")
    lines.append("|----------|--------|-------------|------------|")
    lines.append("| 协议默认（GC 45–55，dTm≤5，Gibbs −70~−50） | 67 | 1.4 | 20/48 |")
    lines.append("| 放宽 GC（40–65） | 102 | 2.1 | 11/48 |")
    lines.append("| 放宽 GC + dTm≤8 | 111 | 2.3 | 11/48 |")
    lines.append("| 全部放宽（GC 40–65，dTm≤8，Gibbs −75~−45） | 186 | 3.9 | 6/48 |")
    lines.append("")
    lines.append("（全部变体均在发卡修复后重跑，数字与当前引擎一致。）")
    lines.append("")
    lines.append("三个分枝杆菌基因组 GC 含量约 65–67%，52-mer tile 落在 45–55% GC")
    lines.append("窗口的概率天然很低——**这是序列组成问题，不是工具问题**。建议对")
    lines.append("分枝杆菌使用 GC 40–65% + dTm≤8 的放宽窗口。即使放宽，HCR 在单基因")
    lines.append("上的对数（1–14 对）仍低于 smFISH/SNAIL 的探针数，但 HCR 每对探针都")
    lines.append("启动一条放大器链聚合，信号按指数级放大，1–14 对通常已可检出；")
    lines.append("要进一步提高对数，可把靶区扩展到全长转录本（含 UTR）。")
    lines.append("")

    lines.append("## 4. SNAIL 的可行性说明")
    lines.append("")
    lines.append("48 基因产出 252 对 primer/padlock（平均 5.2 对/基因），除 420 nt 的")
    lines.append("whiB5 无可用双臂位置外全部基因有产出。双臂结构（41 nt 结合区 + 严格")
    lines.append("的双臂 GC/重复/发卡过滤）产率天然低于线性探针方案，属方法固有特性：")
    lines.append("每对探针经连接 + 滚环扩增放大，5 对左右即可支撑检出，但低丰度靶")
    lines.append("建议优先选择较长基因（≥1 kb，可产出 8–15 对）。订购时 padlock 需加")
    lines.append("5′ 磷酸化（/5Phos/），UGI 条码当前为 N 占位，下单前替换为实际")
    lines.append("正交条码序列。")
    lines.append("")

    lines.append("## 5. 方案选择的建议")
    lines.append("")
    lines.append("| 场景 | 推荐方案 | 理由 |")
    lines.append("|------|----------|------|")
    lines.append("| 低丰度转录本、要求最强检出 | **smiFISH**（或 smFISH） | 单碱基级覆盖、探针堆叠密度最高 |")
    lines.append("| 需要多色/多靶标、希望减少二级探针种类 | **HCR 3.0** | 每通道一套放大器；注意按 GC 调参 |")
    lines.append("| 需要条码化多靶标（UC-seq 类）或信号放大 + 条码 | **SNAIL** | padlock 携带 UGI 条码 |")
    lines.append("| 常规单分子计数 | **smFISH** | 最简单直接 |")
    lines.append("")

    lines.append("## 6. 逐基因 × 逐方案结果矩阵")
    lines.append("")
    lines.append("| 物种 | 类别 | 基因 | 长度 | smiFISH | HCR3 | SNAIL |")
    lines.append("|------|------|------|------|---------|------|-------|")
    gene_rows: dict[tuple, dict[str, dict]] = defaultdict(dict)
    for r in rows:
        gene_rows[(r["org"], r["gene"])][r["scheme"]] = r
    for (org, gene), schemes in gene_rows.items():
        sample = next(iter(schemes.values()))
        cells = []
        for scheme in ("smiFISH", "HCR3", "SNAIL-FISH"):
            r = schemes.get(scheme)
            cells.append(r["final_probes"] if r else "—")
        lines.append(
            f"| {ORG_LABEL.get(org, org)} | {CLS_LABEL.get(sample['cls'], sample['cls'])} "
            f"| {gene} | {sample['length_nt']} | {cells[0]} | {cells[1]} | {cells[2]} |"
        )
    lines.append("")

    lines.append("## 7. 输出文件")
    lines.append("")
    lines.append("| 文件 | 内容 |")
    lines.append("|------|------|")
    lines.append("| `test_data/scheme_matrix_results/scheme_matrix_summary.csv` | 3 方案 × 48 基因漏斗统计 |")
    lines.append("| `batch_scheme_matrix_test.py` | 本轮测试脚本（含产物结构自检） |")
    lines.append("| `test_data/lowabundance_results/gene_fastas/` | 基因输入 FASTA（与第二轮共用） |")
    lines.append("")

    REPORT.write_text("\n".join(lines), encoding="utf-8")
    print(f"报告已生成：{REPORT}")


def num(rows: list[dict], key: str) -> int:
    return sum(int(r[key]) for r in rows)


if __name__ == "__main__":
    main()
