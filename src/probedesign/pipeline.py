"""High-level probe design pipeline."""

from __future__ import annotations

from typing import List

from probedesign.models import DesignParams, DesignResult, ReferenceGenome
from probedesign.schemes import design_for_scheme
from probedesign.schemes.smfish import design_smfish


def run_design(
    target_fasta: str,
    target_index: str,
    host_genomes: List[ReferenceGenome],
    params: DesignParams | None = None,
    threads: int = 1,
) -> DesignResult:
    """Run a probe design job using the scheme selected in ``params``."""
    return design_for_scheme(target_fasta, target_index, host_genomes, params, threads=threads)
