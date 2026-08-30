"""批量探针设计：多记录 FASTA → 逐基因独立设计 → 汇总。

供桌面 GUI 与命令行入口（mycoprimer-batch）共用。

工作方式：
    输入是一个多记录 FASTA（每条记录 = 一个基因/转录本，记录 ID 即基因名）。
    整份输入只构建**一次** bowtie2 索引，随后对每条记录独立跑完整设计流程
    （候选枚举 → 热力学过滤 → 比对 → 打分 → 选点）。共享索引意味着
    target_hits 同时反映"探针在自身基因与其他已提交基因上的额外命中"，
    可用来发现面板内的重复序列。

背景基因组（可选）对所有记录一致：通常**留空**（单菌液培低丰度检测），
仅在需要物种区分时传入（见第一轮测试报告的经验教训）。
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .alignment import build_bowtie2_index
from .models import DesignParams, DesignResult, ReferenceGenome
from .pipeline import run_design

# 项目根目录（CLI 在任意位置运行时用于回退查找注册表）
PROJECT = Path(__file__).resolve().parent.parent.parent
INDICES_DIR = PROJECT / "indices"


@dataclass
class BatchGeneResult:
    """单个基因的批量设计结果。"""

    record_id: str
    description: str
    length_nt: int
    result: DesignResult | None = None  # 失败时为 None
    error: str = ""
    runtime_s: float = 0.0

    @property
    def ok(self) -> bool:
        return self.result is not None

    @property
    def final_probes(self) -> int:
        return len(self.result.passed_probes) if self.result else 0


@dataclass
class BatchReport:
    """整批任务的汇总。"""

    params: DesignParams
    genes: list[BatchGeneResult] = field(default_factory=list)
    elapsed_s: float = 0.0

    @property
    def total_genes(self) -> int:
        return len(self.genes)

    @property
    def total_probes(self) -> int:
        return sum(g.final_probes for g in self.genes)

    @property
    def failed_genes(self) -> list[BatchGeneResult]:
        return [g for g in self.genes if not g.ok]

    def per_gene_rows(self) -> list[dict]:
        """逐基因汇总行（供 CSV 导出）。"""
        rows = []
        for gene in self.genes:
            row = {
                "record_id": gene.record_id,
                "description": gene.description,
                "length_nt": gene.length_nt,
                "final_probes": gene.final_probes,
                "runtime_s": round(gene.runtime_s, 1),
                "error": gene.error,
            }
            if gene.ok and gene.result is not None:
                result = gene.result
                covered = sum(p.stop - p.start for p in result.passed_probes)
                tms = [p.tm for p in result.passed_probes]
                row.update(
                    {
                        "candidates": len(result.probes),
                        "coverage_pct": round(
                            covered / result.target_length * 100, 1
                        ) if result.target_length else 0.0,
                        "tm_mean": round(sum(tms) / len(tms), 1) if tms else "",
                    }
                )
            rows.append(row)
        return rows


def parse_multi_fasta(text: str) -> list[SeqRecord]:
    """解析多记录 FASTA 文本；无记录或存在空序列记录时抛 ValueError。"""
    records = list(SeqIO.parse(io.StringIO(text), "fasta"))
    if not records:
        raise ValueError("FASTA 中没有找到任何序列记录。")
    empty = [r.id or f"记录{i}" for i, r in enumerate(records, 1) if len(r.seq) == 0]
    if empty:
        raise ValueError(f"以下记录没有序列：{', '.join(empty)}。")
    return records


def run_batch(
    fasta_text: str,
    index_prefix: str,
    hosts: list[ReferenceGenome],
    params: DesignParams,
    *,
    cache_dir: Path | None = None,
    threads: int = 2,
    progress: Callable[[int, int, str], None] | None = None,
) -> BatchReport:
    """对多记录 FASTA 的每条记录独立执行一次完整设计。

    参数：
        fasta_text   多记录 FASTA 文本（每条记录 = 一个基因）
        index_prefix 整份输入预构建的 bowtie2 索引前缀（由调用方构建，
                     例如 GUI 的 _build_target_index / CLI 的 build 步骤）
        hosts        背景/宿主基因组（低丰度检测传空列表）
        params       所有记录共用的设计参数
        cache_dir    逐记录临时 FASTA 的写出目录（None 用系统临时目录）
        threads      bowtie2 线程数
        progress     进度回调 progress(当前序号, 总数, 记录 ID)

    返回 BatchReport；单条记录失败不中断整批，错误记入该条 BatchGeneResult。
    """
    records = parse_multi_fasta(fasta_text)
    report = BatchReport(params=params)
    started = time.time()

    tmp_dir = Path(cache_dir) if cache_dir else None
    if tmp_dir is not None:
        tmp_dir.mkdir(parents=True, exist_ok=True)

    for i, record in enumerate(records, start=1):
        if progress is not None:
            progress(i, len(records), record.id)

        fasta_id = record.id or f"record_{i}"
        gene = BatchGeneResult(
            record_id=fasta_id,
            description=record.description,
            length_nt=len(record.seq),
        )
        t0 = time.time()
        if len(record.seq) < params.min_length:
            gene.error = (
                f"记录长度 {len(record.seq)} nt 短于最小探针长度 "
                f"{params.min_length} nt，无法枚举候选。"
            )
            gene.runtime_s = time.time() - t0
            report.genes.append(gene)
            continue
        try:
            record_path = (
                (tmp_dir / f"{_safe_name(fasta_id)}.fasta")
                if tmp_dir
                else Path(_write_temp(fasta_id, record))
            )
            _write_fasta(record_path, record)
            gene.result = run_design(
                str(record_path), index_prefix, hosts, params, threads=threads
            )
        except Exception as exc:  # 单条失败不中断整批
            gene.error = str(exc)[:300]
        gene.runtime_s = time.time() - t0
        report.genes.append(gene)

    report.elapsed_s = time.time() - started
    return report


def _safe_name(record_id: str) -> str:
    keep = "".join(c if c.isalnum() or c in "._-" else "_" for c in record_id)
    return keep or "record"


def _write_temp(record_id: str, record: SeqRecord) -> str:
    import tempfile

    handle = tempfile.NamedTemporaryFile(
        "w", suffix=".fasta", prefix=f"batch_{_safe_name(record_id)}_",
        delete=False, encoding="utf-8",
    )
    _write_fasta(Path(handle.name), record)
    handle.close()
    return handle.name


def _write_fasta(path: Path, record: SeqRecord) -> None:
    seq = str(record.seq).upper()
    with Path(path).open("w", encoding="utf-8") as fh:
        fh.write(f">{record.description or record.id}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i : i + 80] + "\n")


def batch_probe_rows(report: BatchReport) -> list[dict]:
    """全部最终探针合并为行（含所属基因），供 CSV 导出。"""
    rows: list[dict] = []
    for gene in report.genes:
        if not gene.ok or gene.result is None:
            continue
        scheme = gene.result.params.design_scheme
        for probe in gene.result.passed_probes:
            row = {
                "record_id": gene.record_id,
                "probe_id": probe.probe_id,
                "start": probe.start + 1,
                "stop": probe.stop,
                "length": probe.length,
                "sequence": probe.sequence,
                "tm": round(probe.tm, 2),
                "gc": round(probe.gc_content, 3),
                "target_hits": probe.target_hits,
                "score": round(probe.score, 3),
            }
            for key in ("full_sequence", "P1_sequence", "P2_sequence",
                        "primer_sequence", "padlock_sequence_5phos"):
                if key in probe.metadata:
                    row[key] = probe.metadata[key]
            rows.append(row)
    return rows


def main() -> int:
    """命令行批量设计入口：`mycoprimer-batch 输入.fasta [选项]`。"""
    import argparse

    parser = argparse.ArgumentParser(
        prog="mycoprimer-batch",
        description="MycoPrimerV2 批量探针设计：多记录 FASTA 逐基因设计并导出汇总。",
    )
    parser.add_argument("fasta", help="多记录 FASTA 文件路径（每条记录 = 一个基因）")
    parser.add_argument("--scheme", default="smFISH",
                        choices=("smFISH", "smiFISH", "HCR3", "SNAIL-FISH"),
                        help="设计方案（默认 smFISH）")
    parser.add_argument("--desired", type=int, default=48,
                        help="每个基因的目标探针数（默认 48，0 = 不限制）")
    parser.add_argument("--min-gap", type=int, default=2, help="相邻探针最小间隔（默认 2）")
    parser.add_argument("--background", nargs="*", default=[],
                        help="背景基因组 ID 列表（来自注册表；低丰度检测留空）")
    parser.add_argument("--registry", default=None,
                        help="genome_registry.json 路径（默认自动解析）")
    parser.add_argument("--out-dir", default="batch_design_out", help="输出目录")
    parser.add_argument("--threads", type=int, default=2, help="bowtie2 线程数")
    args = parser.parse_args()

    registry_path = (
        Path(args.registry) if args.registry
        else (
            Path(os.environ.get("PROBESTUDIO_HOME", Path.cwd())) / "genome_registry.json"
        )
    )
    if not registry_path.is_file():
        fallback = PROJECT / "genome_registry.json"
        if fallback.is_file():
            registry_path = fallback
    registry = (
        json.loads(registry_path.read_text(encoding="utf-8"))
        if registry_path.is_file() else {}
    )

    hosts = [
        ReferenceGenome(
            id=gid,
            organism=registry.get(gid, {}).get("organism", gid),
            fasta_path=registry.get(gid, {}).get("fasta_path", ""),
            bowtie2_index=registry.get(gid, {}).get("index_prefix", ""),
            is_host=True,
        )
        for gid in args.background if gid in registry
    ]
    unknown = [gid for gid in args.background if gid not in registry]
    if unknown:
        print(f"警告：以下背景基因组未注册，已忽略：{', '.join(unknown)}")

    params = DesignParams(
        design_scheme=args.scheme,
        desired_probe_count=args.desired or None,
        min_gap=args.min_gap,
    )

    fasta_path = Path(args.fasta)
    if not fasta_path.is_file():
        print(f"输入文件不存在：{fasta_path}", file=sys.stderr)
        return 1
    fasta_text = fasta_path.read_text(encoding="utf-8")

    # 整份输入构建一次共享索引（内容哈希缓存）
    INDICES_DIR.mkdir(parents=True, exist_ok=True)
    digest = hashlib.md5(fasta_text.encode()).hexdigest()
    prefix = str(INDICES_DIR / f"_batch_{digest[:12]}")
    if not Path(prefix + ".1.bt2").exists() and not Path(prefix + ".1.bt2l").exists():
        print("正在构建共享 bowtie2 索引…", flush=True)
        build_bowtie2_index(str(fasta_path), prefix, threads=args.threads)

    def progress(i: int, n: int, record_id: str) -> None:
        print(f"  [{i}/{n}] {record_id}", flush=True)

    print(f"批量设计开始：方案 {args.scheme}，"
          f"背景 {', '.join(args.background) or '无'}", flush=True)
    report = run_batch(
        fasta_text, prefix, hosts, params,
        cache_dir=Path(args.out_dir) / "cache",
        threads=args.threads, progress=progress,
    )

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = report.per_gene_rows()
    with (out_dir / "per_gene_summary.csv").open("w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    probe_rows = batch_probe_rows(report)
    with (out_dir / "final_probes.csv").open("w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(probe_rows[0].keys()) if probe_rows else ["record_id"])
        writer.writeheader()
        writer.writerows(probe_rows)

    print(f"完成：{report.total_genes} 个基因 / {report.total_probes} 条探针 / "
          f"{report.elapsed_s:.0f} s；失败 {len(report.failed_genes)} 条。", flush=True)
    print(f"输出目录：{out_dir.resolve()}")
    return 0


if __name__ == "__main__":
    import csv  # noqa: F401  main() 使用
    raise SystemExit(main())
