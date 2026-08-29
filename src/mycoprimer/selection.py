"""筛选后处理：间距筛选与均匀降采样。

select_non_overlapping：
    1. 取出全部通过过滤的候选，按分数从高到低排序；
    2. 贪心遍历：与已选集合在"区间 + min_gap"意义上都不重叠的入选；
       注意是真正的区间重叠判断（v1 只比较起点距离，会放过尾部重叠的探针）；
    3. 若超过 desired_probe_count，再按位置均匀降采样。

equal_space_downsample：
    在候选覆盖的范围 [first.start, last.start] 上等距取 n 个位置，
    各取离目标位置最近的探针，保证最终集合在靶标上均匀铺开。
"""

from __future__ import annotations

from typing import List, Optional

from mycoprimer.models import Probe


def _intervals_too_close(a: Probe, b: Probe, min_gap: int) -> bool:
    """判断两个探针区间是否重叠或间距小于 min_gap（标准区间相交判定）。"""
    return a.start < b.stop + min_gap and b.start < a.stop + min_gap


def select_non_overlapping(
    probes: List[Probe],
    min_gap: int = 0,
    desired_count: Optional[int] = None,
) -> List[Probe]:
    """Greedily select passed probes sorted by score, enforcing min_gap.

    If desired_count is provided and more probes survive spacing than desired,
    evenly downsample the spaced set.
    """
    candidates = [p for p in probes if p.passed]
    candidates.sort(key=lambda p: p.score, reverse=True)

    selected: List[Probe] = []
    for probe in candidates:
        if not any(
            _intervals_too_close(probe, other, min_gap) for other in selected
        ):
            selected.append(probe)

    if desired_count and len(selected) > desired_count:
        selected = equal_space_downsample(selected, desired_count)

    selected.sort(key=lambda p: p.start)
    return selected


def equal_space_downsample(probes: List[Probe], n: int) -> List[Probe]:
    """Pick approximately n evenly spaced probes across the covered region."""
    if n <= 0 or not probes:
        return []
    if len(probes) <= n:
        return probes

    probes = sorted(probes, key=lambda p: p.start)
    if n == 1:
        return [probes[len(probes) // 2]]

    # total span from first start to last start
    first, last = probes[0].start, probes[-1].start
    span = last - first
    if span == 0:
        return probes[:n]

    selected = []
    for i in range(n):
        desired = first + int(round(span * i / (n - 1)))
        closest = min(probes, key=lambda p: abs(p.start - desired))
        if closest not in selected:
            selected.append(closest)

    return sorted(selected, key=lambda p: p.start)
