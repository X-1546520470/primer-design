"""HCR 3.0 probe design scheme."""

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
    """Apply HCR3-specific GC, Tm, hairpin, homopolymer and Gibbs filters."""
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

        hairpin = calc_hairpin(probe.sequence)
        probe.hairpin_tm = hairpin.tm
        if hairpin.tm > params.max_hairpin_tm:
            reasons.append(f"hairpinTm={hairpin.tm:.1f}C > {params.max_hairpin_tm}")

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
    """Split a tile into P1/P2 half-probes and attach HCR initiators."""
    tile = probe.sequence
    tile_len = len(tile)
    mid = tile_len // 2
    five_prime = tile[:mid - 1]
    three_prime = tile[mid + 1:]

    probe.metadata["five_prime_half"] = five_prime
    probe.metadata["three_prime_half"] = three_prime

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
    """Run an HCR 3.0 split-probe design."""
    params = params or DesignParams(design_scheme="HCR3")

    target = maybe_reverse_complement_target(load_first_target(target_fasta), params)
    target_length = len(target.seq)

    # Temporarily override mining length to the HCR tile size.
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
