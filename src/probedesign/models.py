"""Data models for probe design jobs and results."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class ReferenceGenome:
    """一个已注册的参考基因组（用作背景过滤或比对索引）。

    属性说明：
        id             唯一标识，如 "mtb_h37rv"
        organism       展示名称，如 "M. tuberculosis H37Rv"
        fasta_path     基因组 FASTA 文件路径
        bowtie2_index  bowtie2 索引前缀（不含 .1.bt2 等后缀）
        gtf_path       可选的基因注释文件（当前引擎未使用，保留扩展）
        is_host        是否标记为宿主/背景基因组（仅用于界面展示）
    """

    id: str
    organism: str
    fasta_path: str
    bowtie2_index: str
    gtf_path: Optional[str] = None
    is_host: bool = False


@dataclass
class DesignParams:
    """一次探针设计任务的全部参数（所有方案共用一个类，方案只读取自己需要的键）。

    参数分组：
        候选枚举     min_length / max_length / min_tm / max_tm / target_tm
        热力学过滤   min_gc / max_gc / max_homopolymer / max_hairpin_tm
        特异性       bowtie2_preset / bowtie2_score_min / max_target_hits / max_host_hits
        筛选后处理   min_gap（相邻探针最小间隔）/ desired_probe_count（目标条数）
        输入与方案   strand（靶标链向）/ design_scheme（方案名，见 schemes/__init__.py）
        方案专属     smi_*（smiFISH）/ hcr_*（HCR 3.0）/ snail_*（SNAIL FISH）
    """

    # ---- 候选枚举：在靶序列上滑动窗口生成候选 ----
    min_length: int = 18        # 最短窗口（nt）
    max_length: int = 24        # 最长窗口（nt）
    min_tm: float = 50.0        # Tm 下限（°C）
    max_tm: float = 70.0        # Tm 上限（°C）
    target_tm: Optional[float] = 60.0  # 打分时偏好接近该 Tm；None 关闭偏好

    # ---- 热力学过滤 ----
    min_gc: float = 0.20        # GC 分数下限（0-1）
    max_gc: float = 0.80        # GC 分数上限
    max_homopolymer: int = 4    # 允许的最长同聚碱基（超过即过滤，如 GGGGG）
    max_hairpin_tm: float = 45.0  # 发卡 Tm 上限（primer3 计算）

    # ---- 特异性过滤（比对命中数阈值）----
    bowtie2_preset: str = "--very-sensitive-local"  # 灵敏度预设
    bowtie2_score_min: str = "G,20,8"  # 最低比对得分：G,20,8 ≈ 匹配≥20 分
    max_target_hits: int = 10   # 靶标基因组最大命中数（超出视为重复序列）
    max_host_hits: int = 0      # 背景/宿主基因组最大命中数（0 = 零容忍）

    # ---- 筛选后处理 ----
    min_gap: int = 0            # 相邻探针结合区之间的最小间隔（nt）
    desired_probe_count: Optional[int] = None  # 目标探针数（超出则均匀降采样）

    # ---- 输入与方案 ----
    strand: str = "+"           # 靶标链向："+"（mRNA 本身）或 "-"（取反向互补）
    design_scheme: str = "smFISH"  # 设计方案名

    # ---- smiFISH 专属 ----
    smi_readout_sequence: Optional[str] = None  # 共享 FLAP/readout 序列；None 则不拼接
    smi_readout_position: str = "3prime"        # 延伸段位置："3prime" 或 "5prime"
    smi_linker: str = "TTT"                     # 探针与延伸段之间的间隔

    # ---- HCR 3.0 专属 ----
    hcr_tile_size: int = 52      # tile 长度；标准协议 52 nt（两条 25-mer + 中间 2 nt）
    hcr_channel: str = "B1"      # 分裂 initiator 通道（B1–B5），多色实验按靶标区分
    hcr_min_gibbs: float = -70.0  # RNA/DNA 杂交体 Gibbs 自由能窗口下限（kcal/mol）
    hcr_max_gibbs: float = -50.0  # Gibbs 窗口上限
    hcr_dtm_max: Optional[float] = 5.0  # 两条半探针 Tm 差的最大允许值
    hcr_min_gc: float = 45.0     # tile GC 下限（%）
    hcr_max_gc: float = 55.0     # tile GC 上限（%）
    hcr_min_tm: Optional[float] = None  # tile 整体 Tm 窗口下限（默认关闭）
    hcr_max_tm: Optional[float] = None  # tile Tm 上限（52-mer 必然远超 smFISH 窗口）

    # ---- SNAIL FISH 专属 ----
    snail_arm_length: int = 20       # 每条靶结合臂长度（标准 20 nt）
    snail_arm_spacer: int = 1        # 两臂之间在靶标上的间隔 nt 数
    snail_min_gc: float = 40.0       # 单臂 GC 下限（%）
    snail_max_gc: float = 63.0       # 单臂 GC 上限（%）
    snail_hairpin_dg: float = -9.0   # 臂发卡 dG 阈值：低于该值（更负）即过滤
    snail_primer_end: str = "TAATGTTATCTT"   # primer 臂 3′ 端 linker
    snail_padlock_start: str = "ACATTA"      # padlock 5′ anchor
    snail_padlock_end: str = "AAGATA"        # padlock 3′ anchor
    snail_spacer1: str = "ata"       # padlock 内：臂2 与 UGI 条码之间的 spacer
    snail_spacer2: str = "att"       # padlock 内：UGI 条码与 3′ anchor 之间的 spacer
    snail_ugi_sequence: Optional[str] = None  # 正交条码序列；None 则用 N 占位

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
