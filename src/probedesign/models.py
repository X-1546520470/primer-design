"""Data models for probe design jobs and results."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class ReferenceGenome:
    """A registered target or host genome."""

    id: str
    organism: str
    fasta_path: str
    bowtie2_index: str
    gtf_path: Optional[str] = None
    is_host: bool = False


@dataclass
class DesignParams:
    """Parameters for a smFISH probe design run."""

    # Candidate mining
    min_length: int = 18
    max_length: int = 24
    min_tm: float = 50.0
    max_tm: float = 70.0
    target_tm: Optional[float] = 60.0

    # Thermodynamic filters
    min_gc: float = 0.20
    max_gc: float = 0.80
    max_homopolymer: int = 4
    max_hairpin_tm: float = 45.0

    # Specificity
    bowtie2_preset: str = "--very-sensitive-local"
    bowtie2_score_min: str = "G,20,8"
    max_target_hits: int = 10
    max_host_hits: int = 0

    # Selection
    min_gap: int = 0
    desired_probe_count: Optional[int] = None

    # Input handling
    strand: str = "+"
    design_scheme: str = "smFISH"

    # smiFISH
    smi_readout_sequence: Optional[str] = None
    smi_readout_position: str = "3prime"
    smi_linker: str = "TTT"

    # HCR 3.0
    hcr_tile_size: int = 52
    hcr_channel: str = "B1"
    hcr_min_gibbs: float = -70.0
    hcr_max_gibbs: float = -50.0
    hcr_dtm_max: Optional[float] = 5.0
    hcr_min_gc: float = 45.0
    hcr_max_gc: float = 55.0
    hcr_min_tm: Optional[float] = None
    hcr_max_tm: Optional[float] = None

    # SNAIL FISH
    snail_arm_length: int = 20
    snail_arm_spacer: int = 1
    snail_min_gc: float = 40.0
    snail_max_gc: float = 63.0
    snail_hairpin_dg: float = -9.0
    snail_primer_end: str = "TAATGTTATCTT"
    snail_padlock_start: str = "ACATTA"
    snail_padlock_end: str = "AAGATA"
    snail_spacer1: str = "ata"
    snail_spacer2: str = "att"
    snail_ugi_sequence: Optional[str] = None

    def __post_init__(self) -> None:
        if self.min_length > self.max_length:
            raise ValueError("min_length must be <= max_length")
        if self.min_tm > self.max_tm:
            raise ValueError("min_tm must be <= max_tm")
        if self.min_gc > self.max_gc:
            raise ValueError("min_gc must be <= max_gc")


@dataclass
class Probe:
    """A single candidate or final probe."""

    probe_id: str
    target_id: str
    start: int  # 0-based inclusive
    stop: int  # 0-based exclusive
    sequence: str  # target-binding sequence, 5'->3' antisense
    rc_sequence: str  # actual probe sequence

    gc_content: float = 0.0
    tm: float = 0.0
    hairpin_tm: float = 0.0

    target_hits: int = 0
    host_hits: Dict[str, int] = field(default_factory=dict)

    on_target_score: float = 0.0
    off_target_score: float = 0.0
    score: float = 0.0

    passed: bool = True
    failure_reasons: List[str] = field(default_factory=list)
    metadata: Dict = field(default_factory=dict)

    @property
    def length(self) -> int:
        return self.stop - self.start


@dataclass
class DesignResult:
    """Result of a complete design run."""

    params: DesignParams
    target_id: str
    target_length: int
    probes: List[Probe]
    host_genome_ids: List[str]

    @property
    def passed_probes(self) -> List[Probe]:
        return [p for p in self.probes if p.passed]

    @property
    def failed_probes(self) -> List[Probe]:
        return [p for p in self.probes if not p.passed]
