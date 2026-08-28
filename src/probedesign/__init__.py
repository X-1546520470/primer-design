"""ProbeDesign: FISH probe design engine and GUI.

Schemes: smFISH, smiFISH, HCR 3.0, SNAIL-FISH — each in its own module under
``probedesign.schemes``, dispatched by ``design_for_scheme``.
"""

from .models import DesignParams, DesignResult, Probe, ReferenceGenome
from .pipeline import run_design
from .report import probes_to_dataframe, write_outputs

__version__ = "2.0.0"

__all__ = [
    "DesignParams",
    "DesignResult",
    "Probe",
    "ReferenceGenome",
    "probes_to_dataframe",
    "run_design",
    "write_outputs",
]
