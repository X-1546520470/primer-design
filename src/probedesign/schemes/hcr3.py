"""HCR 3.0 探针设计（杂交链式反应放大的 FISH）。

HCR 3.0 的探针结构（Choi et al. 2018 协议）：
    靶标上每 52 nt 的 tile 拆成两条 25-mer 半探针（中间 2 nt 丢弃），
    半探针外侧各带分裂 initiator 的一半：
        P1 = initiator(odd) + 3′ 半探针
        P2 = 5′ 半探针 + initiator(even)
    同通道的放大器发夹识别对应 initiator，启动链式聚合发声。

过滤特点（与 smFISH 不同）：
    - tile GC 窗口 45–55%；
    - RNA/DNA 杂交体 Gibbs 自由能窗口（默认 −70 ~ −50 kcal/mol）；
    - 两条半探针的 dTm ≤ 5 °C（保证两条半链同步杂交）；
    - tile 整体 Tm 不做默认过滤——52-mer 的 Tm 必然远超 smFISH 窗口，
      强行套用会全军覆没（v1 的实际 bug），故 v2 改为可选项。

流程顺序（v2 修正）：热力学过滤 → 拆分与 dTm → 特异性比对 → 打分 → 选点。
拆分必须在选点之前，否则选中的探针可能因 dTm 超标被事后淘汰。
"""

from __future__ import annotations

from typing import List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from primer3 import calc_hairpin, calc_tm

from probedesign.mining import mine_candidates
from probedesign.models import DesignParams, DesignResult, Probe, ReferenceGenome
from probedesign.schemes.common import (
    apply_host_alignment,
    apply_target_alignment,
    calc_gibbs_rna_dna,
    gc_content,
    load_first_target,
    maybe_reverse_complement_target,
)
from probedesign.schemes.initiators import HCR_INITIATORS
from probedesign.scoring import score_probes
from probedesign.selection import select_non_overlapping
from probedesign.utils import calc_tm as probe_calc_tm, has_homopolymer


def _filter_tiles(probes: List[Probe], params: DesignParams) -> None:
    """对 tile 依次应用 GC / Tm(可选) / 同聚碱基 / 发卡 / Gibbs 过滤。"""
    for probe in probes:
        reasons: List[str] = []

        gc = gc_content(probe.sequence)
        probe.gc_content = gc
        if gc < params.hcr_min_gc / 100.0 or gc > params.hcr_max_gc / 100.0:
            reasons.append(
                f"GC={gc:.2f} outside [{params.hcr_min_gc},{params.hcr_max_gc}]"
            )

        tm = probe_calc_tm(probe.sequence)
        probe.tm = tm
        # A 52-mer tile always melts far above smFISH Tm windows; the tile Tm
        # filter is therefore opt-in (the protocol relies on half-probe dTm
        # and the Gibbs window instead).
        if params.hcr_min_tm is not None and params.hcr_max_tm is not None:
            if tm < params.hcr_min_tm or tm > params.hcr_max_tm:
                reasons.append(
                    f"Tm={tm:.1f}C outside [{params.hcr_min_tm},{params.hcr_max_tm}]"
                )

        if has_homopolymer(probe.sequence, params.max_homopolymer):
            reasons.append(f"homopolymer>{params.max_homopolymer}")

        # 注意：发卡检查不在 tile 整体上做——45 °C 阈值是为 18–24-mer 设定的，
        # 52-mer 的发卡 Tm 天然更高（v2 实测会淘汰 ~99% 候选）。半探针层面的
        # 发卡检查移到 _split_and_assemble，针对真正合成的 25-mer 序列。

        gibbs = calc_gibbs_rna_dna(probe.sequence)
        probe.metadata["gibbs_fe"] = gibbs
        if gibbs < params.hcr_min_gibbs or gibbs > params.hcr_max_gibbs:
            reasons.append(
                f"Gibbs={gibbs:.1f} outside [{params.hcr_min_gibbs},{params.hcr_max_gibbs}]"
            )

        if reasons:
            probe.passed = False
            probe.failure_reasons.extend(reasons)


def _split_and_assemble(probe: Probe, params: DesignParams) -> None:
    """把 tile 拆成两条 25-mer 半探针并接上分裂 initiator，同时做 dTm 过滤。

    拆分规则：52-mer 的中间 2 nt 丢弃（协议约定），左右各取 25 nt：
        five_prime  = tile[:25]   → P2 的 5′ 半探针
        three_prime = tile[27:]   → P1 的 3′ 半探针
    dTm = |Tm(five_prime) − Tm(three_prime)|，超过 hcr_dtm_max 则淘汰。
    """
    tile = probe.sequence
    tile_len = len(tile)
    mid = tile_len // 2
    five_prime = tile[:mid - 1]
    three_prime = tile[mid + 1:]

    probe.metadata["five_prime_half"] = five_prime
    probe.metadata["three_prime_half"] = three_prime

    # 发卡检查针对真正合成的 25-mer 半探针（45 °C 阈值在此尺度上才有意义）。
    for label, half in (("5p", five_prime), ("3p", three_prime)):
        hairpin = calc_hairpin(half)
        if hairpin.tm > params.max_hairpin_tm:
            probe.passed = False
            probe.failure_reasons.append(
                f"hairpinTm[{label}]={hairpin.tm:.1f}C > {params.max_hairpin_tm}"
            )
            return

    dtm = abs(calc_tm(five_prime) - calc_tm(three_prime))
    probe.metadata["dTm"] = dtm
    if params.hcr_dtm_max is not None and dtm > params.hcr_dtm_max:
        probe.passed = False
        probe.failure_reasons.append(f"dTm={dtm:.1f} > {params.hcr_dtm_max}")
        return

    initiators = HCR_INITIATORS[params.hcr_channel]
    probe.metadata["channel"] = params.hcr_channel
    probe.metadata["P1_sequence"] = initiators["odd"] + three_prime
    probe.metadata["P2_sequence"] = five_prime + initiators["even"]


def design_hcr3(
    target_fasta: str,
    target_index: str,
    host_genomes: List[ReferenceGenome],
    params: DesignParams | None = None,
    threads: int = 1,
) -> DesignResult:
    """执行一次 HCR 3.0 拆分探针设计。"""
    params = params or DesignParams(design_scheme="HCR3")

    target = maybe_reverse_complement_target(load_first_target(target_fasta), params)
    target_length = len(target.seq)

    # 枚举候选时把窗口长度锁定为 tile 长度（52 nt），其余流程不变。
    from dataclasses import replace
    mining_params = replace(params, min_length=params.hcr_tile_size, max_length=params.hcr_tile_size)
    candidates = mine_candidates(target, mining_params)

    _filter_tiles(candidates, params)

    # Split into halves and apply the dTm filter BEFORE specificity checking
    # and selection, so a selected probe can never be failed afterwards.
    for probe in candidates:
        if probe.passed:
            _split_and_assemble(probe, params)

    apply_target_alignment(candidates, target_index, params, threads)
    apply_host_alignment(candidates, host_genomes, params, threads)

    score_probes(candidates)

    selected = select_non_overlapping(
        candidates,
        min_gap=params.min_gap,
        desired_count=params.desired_probe_count,
    )

    selected_ids = {p.probe_id for p in selected}
    for probe in candidates:
        if probe.passed and probe.probe_id not in selected_ids:
            probe.passed = False
            probe.failure_reasons.append("not_selected")

    return DesignResult(
        params=params,
        target_id=target.id or target.name or "target",
        target_length=target_length,
        probes=candidates,
        host_genome_ids=[h.id for h in host_genomes],
    )
