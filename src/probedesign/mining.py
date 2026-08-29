"""候选探针枚举（mining）。

在靶序列上按 [min_length, max_length] 内的每个长度 L 滑动窗口，枚举全部
(起点, 长度) 组合。例如 1000 nt 靶标、长度 18–24 时约产生 6600 个候选；
随后由 filters / scoring / selection 逐级淘汰。

链向约定：探针的 sequence 字段保存**反义序列**（与靶 RNA 互补配对，
即真正合成使用的序列），因此每个窗口取反向互补后存入。
"""

from __future__ import annotations

from typing import List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from probedesign.models import DesignParams, Probe
from probedesign.utils import reverse_complement


def load_fasta(path: str) -> List[SeqRecord]:
    """读取 FASTA 文件的全部记录（多记录文件只取第一条参与设计）。"""
    return list(SeqIO.parse(path, "fasta"))


def mine_candidates(
    target: SeqRecord,
    params: DesignParams,
) -> List[Probe]:
    """枚举靶序列上的全部候选窗口。

    对每个长度 L ∈ [min_length, max_length]，起点从 0 滑到 len(seq)−L：
        target_window = seq[start:start+L]   # 靶标正链窗口
        probe.sequence = rc(target_window)   # 反义序列 = 合成的探针
        probe.rc_sequence = 正链窗口（备查）

    probe_id 形如 "target:12-36"（0-based，左闭右开）。
    """
    seq = str(target.seq).upper()
    target_id = target.id or target.name or "target"
    candidates: List[Probe] = []

    for length in range(params.min_length, params.max_length + 1):
        for start in range(0, len(seq) - length + 1):
            stop = start + length
            target_window = seq[start:stop]
            probe_seq = reverse_complement(target_window)
            probe_id = f"{target_id}:{start}-{stop}"
            candidates.append(
                Probe(
                    probe_id=probe_id,
                    target_id=target_id,
                    start=start,
                    stop=stop,
                    sequence=probe_seq,
                    rc_sequence=reverse_complement(probe_seq),
                )
            )
    return candidates
