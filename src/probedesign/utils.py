"""Sequence utilities."""

from __future__ import annotations

import math

from probedesign.config import (
    DEFAULT_DNTP,
    DEFAULT_FORMAMIDE_FACTOR,
    DEFAULT_FORMAMIDE_PCT,
    DEFAULT_MG,
    DEFAULT_NA,
    DEFAULT_PROBE_CONC,
    DNA_NN3,
    INITIATION,
    TERMINAL_AT_PENALTY,
)

COMPLEMENT = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")

_GAS_CONSTANT = 1.9872041  # cal / (mol * K)


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA/RNA sequence."""
    return seq.translate(COMPLEMENT)[::-1]


def gc_content(seq: str) -> float:
    """Return GC fraction (0-1) of a sequence."""
    if not seq:
        return 0.0
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq)


def has_homopolymer(seq: str, max_run: int = 4) -> bool:
    """Return True if sequence contains a homopolymer run longer than max_run."""
    if max_run <= 0:
        return False
    seq = seq.upper()
    count = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            count += 1
            if count > max_run:
                return True
        else:
            count = 1
    return False


def salt_correction(
    na: float = DEFAULT_NA, mg: float = DEFAULT_MG, dntp: float = DEFAULT_DNTP
) -> float:
    """Effective monovalent concentration (von Ahsen 2001 Mg2+/dNTP correction)."""
    return na + 120.0 * math.sqrt(max(0.0, mg - dntp))


def calc_tm(
    seq: str,
    probe_conc: float = DEFAULT_PROBE_CONC,
    na: float = DEFAULT_NA,
    mg: float = DEFAULT_MG,
    dntp: float = DEFAULT_DNTP,
    formamide_pct: float = DEFAULT_FORMAMIDE_PCT,
    formamide_factor: float = DEFAULT_FORMAMIDE_FACTOR,
) -> float:
    """Melting temperature via the SantaLucia 1998 nearest-neighbour model.

    Uses the Ct/4 concentration term (non-self-complementary duplex
    assumption) and the linear Mg2+/dNTP salt correction. Formamide lowers
    Tm linearly by ``formamide_factor`` per percent.
    """
    seq = seq.upper()
    if len(seq) == 0:
        return 0.0

    dh, ds = INITIATION

    # terminal AT penalty (positive contribution per SantaLucia/Biopython)
    if seq[0] in ("A", "T"):
        dh += TERMINAL_AT_PENALTY[0]
        ds += TERMINAL_AT_PENALTY[1]
    if seq[-1] in ("A", "T"):
        dh += TERMINAL_AT_PENALTY[0]
        ds += TERMINAL_AT_PENALTY[1]

    for i in range(len(seq) - 1):
        pair = (seq[i], seq[i + 1])
        if pair in DNA_NN3:
            step_dh, step_ds = DNA_NN3[pair]
            dh += step_dh
            ds += step_ds
        else:
            # ambiguous bases: mean stacking values
            dh += -8.0
            ds += -21.0

    monovalent = salt_correction(na, mg, dntp)
    ds += 0.368 * (len(seq) - 1) * math.log(monovalent)

    # Non-self-complementary concentration term: Ct/4.
    tm = (dh * 1000.0) / (ds + _GAS_CONSTANT * math.log(probe_conc / 4.0)) - 273.15
    tm -= formamide_pct * formamide_factor
    return tm


def calc_tm_batch(sequences: dict[str, str], **kwargs) -> dict[str, float]:
    """Calculate Tm for many sequences."""
    return {name: calc_tm(seq, **kwargs) for name, seq in sequences.items()}
