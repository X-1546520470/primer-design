"""Scheme-specific probe design dispatchers."""

from __future__ import annotations

from typing import Dict, List

from probedesign.models import DesignParams, DesignResult, ReferenceGenome
from probedesign.schemes.hcr3 import design_hcr3
from probedesign.schemes.smfish import design_smfish
from probedesign.schemes.smifish import design_smifish
from probedesign.schemes.snail import design_snail

_DESIGN_FUNCS: Dict[str, callable] = {
    "smFISH": design_smfish,
    "smiFISH": design_smifish,
    "HCR3": design_hcr3,
    "SNAIL-FISH": design_snail,
}


def design_for_scheme(
    target_fasta: str,
    target_index: str,
    host_genomes: List[ReferenceGenome],
    params: DesignParams | None = None,
    threads: int = 1,
) -> DesignResult:
    """Dispatch to the appropriate design scheme."""
    params = params or DesignParams()
    scheme = params.design_scheme
    func = _DESIGN_FUNCS.get(scheme)
    if func is None:
        raise ValueError(f"Unknown design scheme: {scheme}")
    return func(target_fasta, target_index, host_genomes, params, threads=threads)
