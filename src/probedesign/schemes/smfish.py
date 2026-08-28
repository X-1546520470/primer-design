"""smFISH probe design scheme."""

from __future__ import annotations

from typing import List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from probedesign.alignment import align_probes_to_index
from probedesign.filters import apply_thermo_filters
from probedesign.mining import load_fasta, mine_candidates
from probedesign.models import DesignParams, DesignResult, Probe, ReferenceGenome
from probedesign.scoring import apply_specificity_filters, score_probes
from probedesign.schemes.common import maybe_reverse_complement_target
from probedesign.selection import select_non_overlapping


def design_smfish(
    target_fasta: str,
    target_index: str,
    host_genomes: List[ReferenceGenome],
    params: DesignParams | None = None,
    threads: int = 1,
) -> DesignResult:
    """Run a complete smFISH probe design with host filtering."""
    params = params or DesignParams(design_scheme="smFISH")

    targets = load_fasta(target_fasta)
    if not targets:
        raise ValueError(f"No sequences found in {target_fasta}")

    target = maybe_reverse_complement_target(targets[0], params)
    target_length = len(target.seq)

    candidates = mine_candidates(target, params)
    apply_thermo_filters(candidates, params)

    target_hits = align_probes_to_index(
        [SeqRecord(seq=Seq(p.sequence), id=p.probe_id, description="") for p in candidates],
        target_index,
        score_min=params.bowtie2_score_min,
        threads=threads,
    )

    host_hits = {}
    for host in host_genomes:
        host_hits[host.id] = align_probes_to_index(
            [SeqRecord(seq=Seq(p.sequence), id=p.probe_id, description="") for p in candidates],
            host.bowtie2_index,
            score_min=params.bowtie2_score_min,
            threads=threads,
        )

    apply_specificity_filters(candidates, target_hits, host_hits, params)
    score_probes(candidates)

    selected = select_non_overlapping(
        candidates,
        min_gap=params.min_gap,
        desired_count=params.desired_probe_count,
    )

    selected_ids = {p.probe_id for p in selected}
    for p in candidates:
        if p.passed and p.probe_id not in selected_ids:
            p.passed = False
            p.failure_reasons.append("not_selected")

    return DesignResult(
        params=params,
        target_id=target.id or target.name or "target",
        target_length=target_length,
        probes=candidates,
        host_genome_ids=[h.id for h in host_genomes],
    )
