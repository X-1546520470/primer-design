"""SNAIL FISH 探针设计（primer + padlock 双臂结构）。

SNAIL（Signal amplification by exchange reaction...）每个靶标位置产出两条寡核苷酸：
    primer  = 臂1（20 nt 靶结合臂）+ 3′ linker
    padlock = 5′ anchor + 臂2（20 nt 靶结合臂）+ spacer + UGI 条码 + spacer + 3′ anchor
两条臂在靶 RNA 上相邻结合（臂间只隔 arm_spacer 个 nt），padlock 的 5′ 磷酸化
末端与 primer 连接后形成环，作为滚环扩增（RCA）模板。

过滤特点：
    - 双臂独立检查 GC（默认 40–63%）/ 重复基序 AAAA·CCCC·GGGG·TTTT /
      发卡 dG（≤ −9 kcal/mol 即淘汰）；
    - 三级特异性检查：整个结合区（cassette）、primer、padlock 分别比对
      靶标与背景基因组，任一级超标即淘汰；
    - 选点时保证相邻探针的双臂互不重叠（最小跨度 = 2×臂长 + 臂间隔）。

链向说明（v2 修正）：结合臂保存为**反义序列**（cassette 的直接切片），
即真正合成并杂交到靶 RNA 的序列。v1 曾把臂还原成正链再组装，得到的
primer/padlock 无法与 RNA 靶标杂交。组装结构本身按项目设计文档保留；
如与你的 SNAIL notebook 有出入，以 notebook 为准并反馈调整。
"""

from __future__ import annotations

from typing import List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from primer3 import calc_hairpin

from probedesign.models import DesignParams, DesignResult, Probe, ReferenceGenome
from probedesign.schemes.common import (
    apply_host_alignment,
    apply_target_alignment,
    calc_hairpin_dg,
    gc_content,
    has_repeat_motif,
    load_first_target,
    maybe_reverse_complement_target,
    reverse_complement,
)
from probedesign.utils import calc_tm
from probedesign.scoring import score_probes
from probedesign.selection import select_non_overlapping


def _build_snail_probe(target_seq: str, start: int, params: DesignParams) -> Probe | None:
    """在 start 位置构造一个 SNAIL 候选。

    靶标正链布局：[臂1(20nt)][间隔(arm_spacer)][臂2(20nt)]
    cassette（探针反义序列）= rc(臂1 + 间隔 + 臂2)，共 41 nt 左右。
    probe.start/stop 覆盖整个双臂区间。
    """
    """Create a SNAIL-FISH candidate with paired arms and full cassette."""
    arm_len = params.snail_arm_length
    spacer = params.snail_arm_spacer
    total = 2 * arm_len + spacer
    if start + total > len(target_seq):
        return None

    arm1 = target_seq[start:start + arm_len]
    arm2 = target_seq[start + arm_len + spacer:start + total]
    cassette = reverse_complement(arm1 + target_seq[start + arm_len:start + arm_len + spacer] + arm2)

    probe_id = f"target:{start}-{start + total}"
    return Probe(
        probe_id=probe_id,
        target_id="target",
        start=start,
        stop=start + total,
        sequence=cassette,
        rc_sequence=reverse_complement(cassette),
    )


def _mine_snail_candidates(target_seq: str, params: DesignParams) -> List[Probe]:
    """按 1 nt 步长枚举全部双臂候选（长度固定为 2×臂长 + 间隔）。"""
    arm_len = params.snail_arm_length
    spacer = params.snail_arm_spacer
    total = 2 * arm_len + spacer
    candidates: List[Probe] = []
    for start in range(0, len(target_seq) - total + 1):
        probe = _build_snail_probe(target_seq, start, params)
        if probe is not None:
            candidates.append(probe)
    return candidates


def _filter_arms(probe: Probe, params: DesignParams) -> bool:
    """Apply per-arm GC, repeat, and hairpin filters. Return True if probe passes.

    The stored arms are the *antisense* (probe-strand) sequences, i.e. direct
    slices of the cassette: these are the sequences actually synthesized and
    the ones that hybridize to the target RNA. The original implementation
    extracted target-strand (sense) arms here, which produced primer/padlock
    oligos unable to bind the RNA target.
    """
    arm_len = params.snail_arm_length
    spacer = params.snail_arm_spacer
    cassette = probe.sequence
    # cassette = rc(arm1_target + spacer + arm2_target), so the antisense arm
    # binding arm1_target is the last arm_len bases and the antisense arm
    # binding arm2_target is the first arm_len bases.
    arm1 = cassette[arm_len + spacer:]
    arm2 = cassette[:arm_len]

    reasons: List[str] = []
    for label, arm in [("arm1", arm1), ("arm2", arm2)]:
        gc = gc_content(arm)
        if gc < params.snail_min_gc / 100.0 or gc > params.snail_max_gc / 100.0:
            reasons.append(f"{label} GC={gc:.2f} outside [{params.snail_min_gc},{params.snail_max_gc}]")
        if has_repeat_motif(arm):
            reasons.append(f"{label} contains AAAA/CCCC/GGGG/TTTT")
        dg = calc_hairpin_dg(arm)
        if dg <= params.snail_hairpin_dg:
            reasons.append(f"{label} hairpin dG={dg:.1f} <= {params.snail_hairpin_dg}")

    if reasons:
        probe.passed = False
        probe.failure_reasons.extend(reasons)
        return False

    probe.metadata["arm1_sequence"] = arm1
    probe.metadata["arm2_sequence"] = arm2
    return True


def _assemble_oligos(probe: Probe, params: DesignParams) -> None:
    """由通过过滤的双臂组装最终订购序列。

    primer  = 臂1(反义) + 3′ linker
    padlock = 5′ anchor + 臂2(反义) + spacer1 + UGI 条码 + spacer2 + 3′ anchor
    UGI 未提供时以 N 占位（订购前需替换为实际条码）；
    padlock 订购时需加 5′ 磷酸化修饰 /5Phos/（GUI 导出表已标注）。
    """
    ugi = params.snail_ugi_sequence or "NNNNNNNNNNNNNNNNNNNNNN"
    primer = probe.metadata["arm1_sequence"] + params.snail_primer_end
    padlock = (
        params.snail_padlock_start
        + probe.metadata["arm2_sequence"]
        + params.snail_spacer1
        + ugi
        + params.snail_spacer2
        + params.snail_padlock_end
    )
    probe.metadata["primer_sequence"] = primer
    probe.metadata["padlock_sequence"] = padlock
    probe.metadata["ugi_barcode"] = ugi


def _check_component_specificity(
    probes: List[Probe],
    sequences: List[str],
    label: str,
    target_index: str,
    host_genomes: List[ReferenceGenome],
    params: DesignParams,
    threads: int,
) -> None:
    """对派生序列（primer / padlock）做独立比对检查（SNAIL 第三级特异性）。

    设计意图：结合区没问题不代表整条 oligo 没问题——linker/anchor/条码
    也可能带来交叉反应。任一组件在背景基因组上超标即淘汰整条候选，
    失败原因记为 "{label}_host_hits[...]" 便于界面区分。
    """
    from probedesign.alignment import align_probes_to_index

    seq_records = [
        SeqRecord(Seq(seq), id=probe.probe_id, description="")
        for probe, seq in zip(probes, sequences)
        if probe.passed
    ]
    if not seq_records:
        return

    target_counts = align_probes_to_index(
        seq_records, target_index, score_min=params.bowtie2_score_min, threads=threads
    )

    host_counts = {}
    for host in host_genomes:
        host_counts[host.id] = align_probes_to_index(
            seq_records, host.bowtie2_index, score_min=params.bowtie2_score_min, threads=threads
        )

    for probe, seq in zip(probes, sequences):
        if not probe.passed:
            continue
        th = target_counts.get(probe.probe_id, 0)
        if th > params.max_target_hits:
            probe.passed = False
            probe.failure_reasons.append(f"{label}_target_hits={th} > {params.max_target_hits}")
            continue
        for host in host_genomes:
            hh = host_counts[host.id].get(probe.probe_id, 0)
            if hh > params.max_host_hits:
                probe.passed = False
                probe.failure_reasons.append(
                    f"{label}_host_hits[{host.id}]={hh} > {params.max_host_hits}"
                )


def design_snail(
    target_fasta: str,
    target_index: str,
    host_genomes: List[ReferenceGenome],
    params: DesignParams | None = None,
    threads: int = 1,
) -> DesignResult:
    """执行一次 SNAIL FISH 设计（双臂过滤 + 三级特异性 + 双臂间距选点）。"""
    params = params or DesignParams(design_scheme="SNAIL-FISH")

    target = maybe_reverse_complement_target(load_first_target(target_fasta), params)
    target_seq = str(target.seq).upper()
    target_length = len(target_seq)

    candidates = _mine_snail_candidates(target_seq, params)
    for probe in candidates:
        probe.tm = calc_tm(probe.sequence)
        probe.gc_content = gc_content(probe.sequence)
        _filter_arms(probe, params)

    # Three-level specificity: full cassette, primer, padlock.
    apply_target_alignment(candidates, target_index, params, threads)
    apply_host_alignment(candidates, host_genomes, params, threads)

    # Assemble oligos for still-passing candidates before checking them.
    for probe in candidates:
        if probe.passed:
            _assemble_oligos(probe, params)

    primer_seqs = [p.metadata.get("primer_sequence", "") for p in candidates]
    padlock_seqs = [p.metadata.get("padlock_sequence", "") for p in candidates]
    _check_component_specificity(
        candidates, primer_seqs, "primer", target_index, host_genomes, params, threads
    )
    _check_component_specificity(
        candidates, padlock_seqs, "padlock", target_index, host_genomes, params, threads
    )

    score_probes(candidates)

    # SNAIL spacing: adjacent arm pairs must not overlap.
    arm_len = params.snail_arm_length
    spacer = params.snail_arm_spacer
    min_span = 2 * arm_len + spacer + params.min_gap
    selected = select_non_overlapping(
        candidates,
        min_gap=min_span,
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
