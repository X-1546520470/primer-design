"""热力学与序列复杂度过滤。

过滤顺序（性能关键）：
    1. GC 分数        —— 纯字符统计，近乎零开销
    2. 同聚碱基       —— 线性扫描，开销极低
    3. Tm             —— 最近邻累加，开销小
    4. 发卡折叠       —— primer3 热力学计算，单条 ~0.1–1 ms，最贵

因此前三级全部通过的候选才送去折叠计算；在数千候选规模下可比"一次性全算"
快一个数量级。未通过前三级的原因也会记入 failure_reasons，便于界面展示。
"""

from __future__ import annotations

from typing import List

from primer3 import calc_hairpin

from probedesign.models import DesignParams, Probe
from probedesign.utils import calc_tm, gc_content, has_homopolymer


def apply_thermo_filters(probes: List[Probe], params: DesignParams) -> List[Probe]:
    """对候选列表依次应用 GC / 同聚碱基 / Tm / 发卡过滤。

    直接修改每个 probe 的 passed 与 failure_reasons 字段（就地修改），
    同时把测得的 gc_content / tm / hairpin_tm 写回探针供展示与打分。
    返回同一列表以便链式调用。
    """
    for probe in probes:
        reasons: List[str] = []

        # Cheap sequence-level checks first; the expensive primer3 hairpin
        # folding runs only on candidates that survive them.
        gc = gc_content(probe.sequence)
        probe.gc_content = gc
        if gc < params.min_gc or gc > params.max_gc:
            reasons.append(f"GC={gc:.2f} outside [{params.min_gc:.2f},{params.max_gc:.2f}]")

        if has_homopolymer(probe.sequence, params.max_homopolymer):
            reasons.append(f"homopolymer>{params.max_homopolymer}")

        tm = calc_tm(probe.sequence)
        probe.tm = tm
        if tm < params.min_tm or tm > params.max_tm:
            reasons.append(f"Tm={tm:.1f}C outside [{params.min_tm:.1f},{params.max_tm:.1f}]")

        if reasons:
            probe.passed = False
            probe.failure_reasons.extend(reasons)
            continue

        hairpin = calc_hairpin(probe.sequence)
        probe.hairpin_tm = hairpin.tm
        if hairpin.tm > params.max_hairpin_tm:
            probe.passed = False
            probe.failure_reasons.append(
                f"hairpinTm={hairpin.tm:.1f}C > {params.max_hairpin_tm:.1f}"
            )

    return probes
