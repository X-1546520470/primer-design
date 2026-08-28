"""Thermodynamic and sequence-complexity filters."""

from __future__ import annotations

from typing import List

from primer3 import calc_hairpin

from probedesign.models import DesignParams, Probe
from probedesign.utils import calc_tm, gc_content, has_homopolymer


def apply_thermo_filters(probes: List[Probe], params: DesignParams) -> List[Probe]:
    """Apply GC, Tm, homopolymer and hairpin filters to candidates.

    Mutates probe.passed and probe.failure_reasons in place.
    Returns the same list for convenience.
    """
    for probe in probes:
        reasons: List[str] = []

        # Cheap sequence-level checks first; the expensive primer3 hairpin
        # folding runs only on candidates that survive them.
        gc = gc_content(probe.sequence)
        probe.gc_content = gc
        if gc < params.min_gc or gc > params.max_gc:
            reasons.append(f"GC={gc:.2f} outside [{params.min_gc:.2f},{params.max_gc:.2f}]")

        if has_homopolymer(probe.sequence, params.max_homopolymer):
            reasons.append(f"homopolymer>{params.max_homopolymer}")

        tm = calc_tm(probe.sequence)
        probe.tm = tm
        if tm < params.min_tm or tm > params.max_tm:
            reasons.append(f"Tm={tm:.1f}C outside [{params.min_tm:.1f},{params.max_tm:.1f}]")

        if reasons:
            probe.passed = False
            probe.failure_reasons.extend(reasons)
            continue

        hairpin = calc_hairpin(probe.sequence)
        probe.hairpin_tm = hairpin.tm
        if hairpin.tm > params.max_hairpin_tm:
            probe.passed = False
            probe.failure_reasons.append(
                f"hairpinTm={hairpin.tm:.1f}C > {params.max_hairpin_tm:.1f}"
            )

    return probes
