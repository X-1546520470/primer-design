"""Candidate probe mining from target sequences."""

from __future__ import annotations

from typing import List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from probedesign.models import DesignParams, Probe
from probedesign.utils import reverse_complement


def load_fasta(path: str) -> List[SeqRecord]:
    """Load all records from a FASTA file."""
    return list(SeqIO.parse(path, "fasta"))


def mine_candidates(
    target: SeqRecord,
    params: DesignParams,
) -> List[Probe]:
    """Generate all candidate probes from a target sequence.

    Probes are extracted as reverse-complement windows so that the
    returned `sequence` is antisense to the target (i.e., the probe
    sequence that will hybridize to the target RNA/DNA).
    """
    seq = str(target.seq).upper()
    target_id = target.id or target.name or "target"
    candidates: List[Probe] = []

    for length in range(params.min_length, params.max_length + 1):
        for start in range(0, len(seq) - length + 1):
            stop = start + length
            target_window = seq[start:stop]
            probe_seq = reverse_complement(target_window)
            probe_id = f"{target_id}:{start}-{stop}"
            candidates.append(
                Probe(
                    probe_id=probe_id,
                    target_id=target_id,
                    start=start,
                    stop=stop,
                    sequence=probe_seq,
                    rc_sequence=reverse_complement(probe_seq),
                )
            )
    return candidates
