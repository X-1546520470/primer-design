"""Default configurations and constants."""

from __future__ import annotations

from probedesign.models import DesignParams

DEFAULT_PARAMS = DesignParams()

# SantaLucia 1998 unified nearest-neighbour table (DNA_NN3), indexed by the
# 5'->3' dinucleotide of the probe strand. dH in kcal/mol, dS in cal/(mol*K).
# All 16 dinucleotides present; symmetric steps (AA/TT, GG/CC) share values.
DNA_NN3 = {
    ("A", "A"): (-7.9, -22.2),
    ("T", "T"): (-7.9, -22.2),
    ("A", "T"): (-7.2, -20.4),
    ("T", "A"): (-7.2, -21.3),
    ("A", "G"): (-7.8, -21.0),  # AG/CT
    ("C", "T"): (-7.8, -21.0),
    ("G", "A"): (-8.2, -22.2),  # GA/CT
    ("T", "C"): (-8.2, -22.2),
    ("G", "G"): (-8.0, -19.9),
    ("C", "C"): (-8.0, -19.9),
    ("G", "C"): (-10.6, -27.2),
    ("C", "G"): (-9.8, -24.4),
    ("A", "C"): (-8.4, -22.4),  # AC/GT
    ("G", "T"): (-8.4, -22.4),
    ("C", "A"): (-8.5, -22.7),  # CA/GT
    ("T", "G"): (-8.5, -22.7),
}

# Terminal A/T initiation penalty (Biopython DNA_NN3 convention: +dH, +dS
# per terminal A/T base pair). Applied to the 5' and 3' terminal bases.
TERMINAL_AT_PENALTY = (2.2, 6.9)

# Duplex initiation values (SantaLucia 1998).
INITIATION = (0.2, -5.7)

# Default salt / formamide conditions for Tm. 0.39 M monovalent approximates
# 2x SSC, the standard FISH hybridization buffer.
DEFAULT_NA = 0.39  # M
DEFAULT_MG = 0.0  # M
DEFAULT_DNTP = 0.0  # M
DEFAULT_PROBE_CONC = 1e-6  # M
DEFAULT_FORMAMIDE_PCT = 0.0
DEFAULT_FORMAMIDE_FACTOR = 0.65  # degC per percent formamide
