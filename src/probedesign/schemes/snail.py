"""SNAIL-FISH probe design scheme.

SNAIL-FISH produces a primer oligo and a 5'-phosphorylated padlock oligo for
each target position. The two 20-nt target-binding arms are separated by a
short spacer on the target RNA.
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
    """Mine all possible SNAIL-FISH paired-arm candidates."""
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
    """Build full primer and padlock oligos from a passing SNAIL candidate."""
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
    """Align a set of derived sequences (primer/padlock/cassette) and mark failures."""
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
    """Run a SNAIL-FISH primer + padlock design."""
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
