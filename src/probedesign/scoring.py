"""Probe scoring and specificity filtering."""

from __future__ import annotations

from typing import Dict, List

from probedesign.models import DesignParams, Probe


def apply_specificity_filters(
    probes: List[Probe],
    target_hits: Dict[str, int],
    host_hits: Dict[str, Dict[str, int]],
    params: DesignParams,
) -> List[Probe]:
    """Apply target-genome and host-genome hit filters.

    Mutates probes in place. Returns the same list.
    """
    for probe in probes:
        if not probe.passed:
            continue

        hits = target_hits.get(probe.probe_id, 0)
        probe.target_hits = hits
        if hits > params.max_target_hits:
            probe.passed = False
            probe.failure_reasons.append(f"target_hits={hits} > {params.max_target_hits}")
            continue

        for host_id, counts in host_hits.items():
            host_h = counts.get(probe.probe_id, 0)
            probe.host_hits[host_id] = host_h
            if host_h > params.max_host_hits:
                probe.passed = False
                probe.failure_reasons.append(f"host_hits[{host_id}]={host_h} > {params.max_host_hits}")

    return probes


def score_probes(probes: List[Probe]) -> List[Probe]:
    """Compute a simple score for ranking.

    Higher is better. Rewards low off-target hits and Tm near 60 C.
    """
    for probe in probes:
        total_host_hits = sum(probe.host_hits.values())
        probe.off_target_score = probe.target_hits + total_host_hits + 1
        probe.on_target_score = 1.0
        # Prefer Tm near 60; mild penalty for extreme GC
        tm_penalty = abs(probe.tm - 60.0) / 20.0
        gc_penalty = abs(probe.gc_content - 0.5) * 2.0
        probe.score = (probe.on_target_score / probe.off_target_score) - tm_penalty - gc_penalty
    return probes
