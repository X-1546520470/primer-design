"""各设计方案共享的辅助函数。

内容分三组：
    1. RNA/DNA 杂交体 Gibbs 自由能（Sugimoto 1995 参数）—— HCR 3.0 的
       核心过滤指标；盐修正沿用 SantaLucia 式线性项；
    2. 靶序列加载与链向处理（strand == '-' 时取反向互补）；
    3. 比对过滤封装：把候选（或 SNAIL 派生的 primer/padlock）比对到
       靶标与背景基因组并按阈值淘汰。
"""

from __future__ import annotations

import math
from typing import List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from primer3 import calc_hairpin

from probedesign.models import DesignParams, Probe, ReferenceGenome
from probedesign.utils import gc_content, reverse_complement


# Sugimoto 1995 RNA/DNA 杂交体最近邻参数（键为小写二核苷酸）。
# dH 单位 kcal/mol，dS 单位 cal/(mol*K)。用于估算探针与靶 RNA 杂交体的
# 结合自由能（HCR 官方推荐窗口约 −50 ~ −70 kcal/mol）。
_SUGIMOTO_DEL_H = {
    "aa": -7.8, "ac": -5.9, "ag": -9.1, "at": -8.3,
    "ca": -9.0, "cc": -9.3, "cg": -16.3, "ct": -7.0,
    "ga": -5.5, "gc": -8.0, "gg": -12.8, "gt": -7.8,
    "ta": -7.8, "tc": -8.6, "tg": -10.4, "tt": -11.5,
}
_SUGIMOTO_DEL_S = {
    "aa": -21.9, "ac": -12.3, "ag": -23.5, "at": -23.9,
    "ca": -26.1, "cc": -23.2, "cg": -47.1, "ct": -19.7,
    "ga": -13.5, "gc": -17.1, "gg": -31.9, "gt": -21.6,
    "ta": -23.2, "tc": -22.9, "tg": -28.4, "tt": -36.4,
}


def calc_gibbs_rna_dna(sequence: str, temp_c: float = 37.0, salt_m: float = 0.33) -> float:
    """计算 RNA/DNA 杂交体的 Gibbs 自由能（kcal/mol），值越负结合越强。

        ΔG = ΔH − T·ΔS + 盐修正(−0.114·N·ln[Na⁺])

    参数 sequence 为 DNA 探针序列（与 RNA 靶互补）；杂交体按 RNA/DNA
    双链处理，采用 Sugimoto 1995 最近邻参数。HCR 3.0 用它筛选
    "结合强度适中"的 tile：太弱杂交不上，太强则背景高。
    """
    seq = sequence.lower()
    if len(seq) < 2:
        return 0.0

    dH = sum(_SUGIMOTO_DEL_H[seq[i:i + 2]] for i in range(len(seq) - 1))
    dS = sum(_SUGIMOTO_DEL_S[seq[i:i + 2]] for i in range(len(seq) - 1))

    # Initiation values from Sugimoto et al. 1995.
    dH += 1.9
    dS += -3.9

    # Gibbs in cal/mol, then convert to kcal/mol.
    g_cal = dH * 1000.0 - (temp_c + 273.15) * dS
    g_kcal = g_cal / 1000.0

    # Salt correction (SantaLucia 1998, eq. 7).
    return g_kcal - 0.114 * len(seq) * math.log(salt_m)


def calc_hairpin_dg(sequence: str) -> float:
    """primer3 计算的发卡自由能（kcal/mol）；SNAIL 双臂过滤用。"""
    result = calc_hairpin(sequence)
    return result.dg / 1000.0 if hasattr(result, "dg") else 0.0


def has_repeat_motif(sequence: str, motifs: List[str] | None = None) -> bool:
    """检测是否含有禁用重复基序；SNAIL 标准为 AAAA/CCCC/GGGG/TTTT。"""
    if motifs is None:
        motifs = ["AAAA", "CCCC", "GGGG", "TTTT"]
    seq = sequence.upper()
    return any(motif in seq for motif in motifs)


def load_first_target(target_fasta: str) -> SeqRecord:
    """Load the first record from a FASTA file."""
    from Bio import SeqIO
    records = list(SeqIO.parse(target_fasta, "fasta"))
    if not records:
        raise ValueError(f"No sequences found in {target_fasta}")
    return records[0]


def maybe_reverse_complement_target(target: SeqRecord, params: DesignParams) -> SeqRecord:
    """Return the target, reverse-complemented if params.strand == '-'."""
    if params.strand == "-":
        rc = target.reverse_complement()
        rc.id = target.id or target.name or "target"
        return rc
    return target


def apply_host_alignment(
    probes: List[Probe],
    host_genomes: List[ReferenceGenome],
    params: DesignParams,
    threads: int,
) -> None:
    """Align probes to each host genome and record host hits.

    Mutates probes in place. Skips probes that already failed.
    """
    from probedesign.alignment import align_probes_to_index

    seq_records = [
        SeqRecord(Seq(p.sequence), id=p.probe_id, description="") for p in probes if p.passed
    ]
    if not seq_records:
        return

    for host in host_genomes:
        counts = align_probes_to_index(
            seq_records, host.bowtie2_index, score_min=params.bowtie2_score_min, threads=threads
        )
        for probe in probes:
            if not probe.passed:
                continue
            host_h = counts.get(probe.probe_id, 0)
            probe.host_hits[host.id] = host_h
            if host_h > params.max_host_hits:
                probe.passed = False
                probe.failure_reasons.append(f"host_hits[{host.id}]={host_h} > {params.max_host_hits}")


def apply_target_alignment(
    probes: List[Probe],
    target_index: str,
    params: DesignParams,
    threads: int,
) -> None:
    """Align probes to target genome and record target hits, filtering by max_target_hits."""
    from probedesign.alignment import align_probes_to_index

    seq_records = [
        SeqRecord(Seq(p.sequence), id=p.probe_id, description="") for p in probes if p.passed
    ]
    if not seq_records:
        return

    counts = align_probes_to_index(
        seq_records, target_index, score_min=params.bowtie2_score_min, threads=threads
    )
    for probe in probes:
        if not probe.passed:
            continue
        hits = counts.get(probe.probe_id, 0)
        probe.target_hits = hits
        if hits > params.max_target_hits:
            probe.passed = False
            probe.failure_reasons.append(f"target_hits={hits} > {params.max_target_hits}")
