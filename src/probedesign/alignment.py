"""Bowtie2 alignment wrappers and SAM parsing."""

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
    """Align sequences to a Bowtie2 index and return hit counts per sequence ID.

    Parameters
    ----------
    sequences : List[SeqRecord]
        Sequences to align.
    index_prefix : str
        Path prefix of the Bowtie2 index.
    score_min : str
        Bowtie2 --score-min argument.
    preset : str
        Bowtie2 sensitivity preset.
    k : int
        Report up to k alignments per read.
    threads : int
        Bowtie2 threads.

    Returns
    -------
    Dict[str, int]
        Mapping from sequence ID to number of reported alignments.
        Unmapped sequences get 0.
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
    """Convenience wrapper that returns hit counts indexed by probe_id."""
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
