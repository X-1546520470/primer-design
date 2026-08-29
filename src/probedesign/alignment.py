"""bowtie2 比对封装与 SAM 命中计数。

特异性过滤的基础设施：把候选探针比对到靶标/背景基因组的 bowtie2 索引，
统计每条探针的命中数。

关键实现细节：
    - 比对参数：--very-sensitive-local + --score-min G,20,8（≈ 匹配 ≥20 分，
      即 18–24 nt 探针容忍少量错配）+ -k 100（最多报告 100 个比对位点）；
    - 命中计数必须**包含 secondary 比对**（SAM flag 0x100）：-k 模式下每个
      位点一条记录，只有一条是 primary，其余全部带 0x100 标志。v1 的解析器
      把 secondary 全部跳过，导致每条探针最多计 1 次命中，重复序列与宿主
      过滤实际失效——这是 v2 修复的最重要 bug；
    - 探针 ID 含冒号，bowtie2 的 QNAME 处理不可靠，因此提交比对前先把
      ID 映射为纯数字，比对完再映射回来；
    - bowtie2 可执行文件优先从当前 conda 环境的 bin 目录查找（sys.executable
      同级），找不到再查 PATH；全部错误包装为 AlignmentError 并保留 stderr。
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class AlignmentError(RuntimeError):
    """Raised when bowtie2 is missing or an alignment run fails."""


def _find_binary(name: str) -> str:
    """Locate a bowtie2-family binary: env bin dir first, then PATH."""
    candidate = Path(sys.executable).resolve().parent / name
    if candidate.is_file() and os.access(candidate, os.X_OK):
        return str(candidate)
    discovered = shutil.which(name)
    if discovered:
        return discovered
    raise AlignmentError(
        f"找不到 {name}。请确认已安装（conda install -c bioconda bowtie2）"
        "且位于当前 conda 环境或 PATH 中。"
    )


def run_bowtie2(
    sequences: List[SeqRecord],
    index_prefix: str,
    score_min: str = "G,20,8",
    preset: str = "--very-sensitive-local",
    k: int = 100,
    threads: int = 1,
) -> Dict[str, int]:
    """把序列比对到 bowtie2 索引，返回每条序列的命中数。

    参数：
        sequences    待比对的序列记录（ID 会映射为数字后提交）
        index_prefix bowtie2 索引前缀（不含 .1.bt2 等后缀）
        score_min    --score-min 参数，默认 G,20,8
        preset       灵敏度预设，默认 --very-sensitive-local
        k            每条读取最多报告的比对位点数
        threads      bowtie2 线程数

    返回：{序列 ID: 命中位点数}；未比对上的为 0。
    比对在临时目录进行，FASTA 与 SAM 用完即删。
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, "queries.fa")
        sam_path = os.path.join(tmpdir, "alignments.sam")
        SeqIO.write(sequences, fasta_path, "fasta")

        cmd = [
            _find_binary("bowtie2"),
            preset,
            "-f",
            "--no-sq",
            "--no-hd",
            "--reorder",
            "--score-min",
            score_min,
            "-k",
            str(k),
            "-p",
            str(threads),
            "-x",
            index_prefix,
            "-U",
            fasta_path,
            "-S",
            sam_path,
        ]
        try:
            completed = subprocess.run(
                cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True
            )
        except FileNotFoundError as exc:
            raise AlignmentError(f"无法启动 bowtie2：{exc}") from exc
        except subprocess.CalledProcessError as exc:
            detail = (exc.stderr or "").strip() or "无错误详情"
            raise AlignmentError(f"bowtie2 运行失败：{detail}") from exc

        return parse_sam_hit_counts(sam_path, expected=len(sequences))


def parse_sam_hit_counts(sam_path: str, expected: int = 0) -> Dict[str, int]:
    """Parse a simple SAM file and count alignments per query ID.

    Counts *all* reported alignments per query, including the secondary
    (flag 0x100) records emitted by ``-k`` mode: bowtie2 reports one primary
    plus up to k-1 secondary alignments, and the total is the per-probe hit
    count used for specificity filtering. Unmapped queries keep 0.
    """
    counts: Dict[str, int] = {str(i): 0 for i in range(1, expected + 1)}
    with open(sam_path, "r") as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 11:
                continue
            qname, flag, rname = parts[0], int(parts[1]), parts[2]
            if rname == "*":
                continue
            if flag & 0x4:  # segment unmapped
                continue
            counts[qname] = counts.get(qname, 0) + 1
    return counts


def build_bowtie2_index(fasta_path: str, index_prefix: str, threads: int = 1) -> None:
    """Build a Bowtie2 index from a FASTA file.

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file.
    index_prefix : str
        Output index prefix.
    threads : int
        Threads for bowtie2-build.
    """
    if not os.path.isfile(fasta_path):
        raise AlignmentError(f"FASTA 文件不存在：{fasta_path}")
    cmd = [
        _find_binary("bowtie2-build"),
        "--threads",
        str(threads),
        fasta_path,
        index_prefix,
    ]
    try:
        subprocess.run(
            cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True
        )
    except FileNotFoundError as exc:
        raise AlignmentError(f"无法启动 bowtie2-build：{exc}") from exc
    except subprocess.CalledProcessError as exc:
        detail = (exc.stderr or "").strip() or "无错误详情"
        raise AlignmentError(f"bowtie2-build 运行失败：{detail}") from exc


def align_probes_to_index(
    probes: List[SeqRecord],
    index_prefix: str,
    score_min: str = "G,20,8",
    threads: int = 1,
) -> Dict[str, int]:
    """便捷封装：提交探针列表并返回 {probe_id: 命中数}，索引缺失时给出明确报错。"""
    if not os.path.exists(index_prefix + ".1.bt2") and not os.path.exists(
        index_prefix + ".1.bt2l"
    ):
        raise AlignmentError(
            f"Bowtie2 索引不存在：{index_prefix}。请先构建索引。"
        )
    # Bowtie2 qnames must not contain colons in some contexts; use a numeric map.
    indexed: List[SeqRecord] = []
    id_map: Dict[str, str] = {}
    for i, rec in enumerate(probes, 1):
        numeric_id = str(i)
        id_map[numeric_id] = rec.id
        indexed.append(SeqRecord(Seq(str(rec.seq)), id=numeric_id, description=""))

    raw_counts = run_bowtie2(
        indexed,
        index_prefix,
        score_min=score_min,
        threads=threads,
        k=100,
    )
    return {id_map[qid]: raw_counts.get(qid, 0) for qid in id_map}
