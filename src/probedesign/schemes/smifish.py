"""smiFISH probe design scheme.

smiFISH is smFISH plus a shared terminal readout extension that binds a
fluorescent secondary probe.
"""

from __future__ import annotations

from typing import List

from probedesign.models import DesignParams, DesignResult, ReferenceGenome
from probedesign.schemes.smfish import design_smfish


def design_smifish(
    target_fasta: str,
    target_index: str,
    host_genomes: List[ReferenceGenome],
    params: DesignParams | None = None,
    threads: int = 1,
) -> DesignResult:
    """Run smFISH design, then append a shared readout extension to all probes."""
    params = params or DesignParams(design_scheme="smiFISH")
    result = design_smfish(target_fasta, target_index, host_genomes, params, threads=threads)

    readout = (params.smi_readout_sequence or "").upper()
    linker = params.smi_linker.upper()
    position = params.smi_readout_position

    for probe in result.probes:
        full_seq = probe.sequence
        if readout:
            if position == "5prime":
                full_seq = readout + linker + probe.sequence
            else:
                full_seq = probe.sequence + linker + readout
        probe.metadata["full_sequence"] = full_seq
        probe.metadata["readout_sequence"] = readout
        probe.metadata["readout_position"] = position
        probe.metadata["linker"] = linker

    return result
