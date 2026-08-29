"""smFISH 探针设计（最基础的方案，也是 smiFISH 的底层）。

设计流程 design_smfish：
    1. 读入靶序列（按 strand 决定是否反向互补）；
    2. mining：滑动窗口枚举候选（长度 min_length–max_length）；
    3. filters：GC / 同聚碱基 / Tm / 发卡热力学过滤；
    4. 比对：候选比对到靶标基因组（检测重复）与各背景基因组（宿主过滤）；
    5. scoring：按特异性 + Tm/GC 偏好打分；
    6. selection：按分数贪心选取互不重叠的最终集合（可限数量）。

产物：每条探针一条反义寡核苷酸，sequence 字段即为订购序列。
"""

from __future__ import annotations

from typing import List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from mycoprimer.alignment import align_probes_to_index
from mycoprimer.filters import apply_thermo_filters
from mycoprimer.mining import load_fasta, mine_candidates
from mycoprimer.models import DesignParams, DesignResult, Probe, ReferenceGenome
from mycoprimer.scoring import apply_specificity_filters, score_probes
from mycoprimer.schemes.common import maybe_reverse_complement_target
from mycoprimer.selection import select_non_overlapping


def design_smfish(
    target_fasta: str,
    target_index: str,
    host_genomes: List[ReferenceGenome],
    params: DesignParams | None = None,
    threads: int = 1,
) -> DesignResult:
    """执行一次完整的 smFISH 设计（含宿主基因组过滤）。"""
    params = params or DesignParams(design_scheme="smFISH")

    targets = load_fasta(target_fasta)
    if not targets:
        raise ValueError(f"No sequences found in {target_fasta}")

    target = maybe_reverse_complement_target(targets[0], params)
    target_length = len(target.seq)

    candidates = mine_candidates(target, params)
    apply_thermo_filters(candidates, params)

    target_hits = align_probes_to_index(
        [SeqRecord(seq=Seq(p.sequence), id=p.probe_id, description="") for p in candidates],
        target_index,
        score_min=params.bowtie2_score_min,
        threads=threads,
    )

    host_hits = {}
    for host in host_genomes:
        host_hits[host.id] = align_probes_to_index(
            [SeqRecord(seq=Seq(p.sequence), id=p.probe_id, description="") for p in candidates],
            host.bowtie2_index,
            score_min=params.bowtie2_score_min,
            threads=threads,
        )

    apply_specificity_filters(candidates, target_hits, host_hits, params)
    score_probes(candidates)

    selected = select_non_overlapping(
        candidates,
        min_gap=params.min_gap,
        desired_count=params.desired_probe_count,
    )

    selected_ids = {p.probe_id for p in selected}
    for p in candidates:
        if p.passed and p.probe_id not in selected_ids:
            p.passed = False
            p.failure_reasons.append("not_selected")

    return DesignResult(
        params=params,
        target_id=target.id or target.name or "target",
        target_length=target_length,
        probes=candidates,
        host_genome_ids=[h.id for h in host_genomes],
    )
