"""批量设计模块（mycoprimer.batch）的单元测试。"""

from __future__ import annotations

import unittest
from pathlib import Path

from mycoprimer.batch import (
    BatchReport,
    batch_probe_rows,
    parse_multi_fasta,
    run_batch,
)
from mycoprimer.models import DesignParams, ReferenceGenome

PROJECT = Path(__file__).resolve().parent.parent
SMOKE_INDEX = PROJECT / "data" / "smoke" / "target_idx"

GENE_A = "ATGACCATGATTACGCCAAGCGCGCTTTTTGCGCGCGATTACAGATTACAGATTAC"
GENE_B = "TTGACCGATGACCCCGGTTCAGGCTTCACCACAGTGTGGAACGCGGTCGTCTCCGAACTT"


def _two_gene_fasta() -> str:
    return (
        ">geneA synthetic\n" + GENE_A + "\n"
        ">geneB synthetic\n" + GENE_B + "\n"
    )


class ParseMultiFastaTests(unittest.TestCase):
    def test_parses_two_records(self) -> None:
        records = parse_multi_fasta(_two_gene_fasta())
        self.assertEqual([r.id for r in records], ["geneA", "geneB"])
        self.assertEqual(len(records[0].seq), len(GENE_A))

    def test_empty_input_raises(self) -> None:
        with self.assertRaises(ValueError):
            parse_multi_fasta(">only_header_no_seq\n")


class RunBatchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not SMOKE_INDEX.with_suffix("").exists() and not (
            PROJECT / "data/smoke/target_idx.1.bt2"
        ).exists():
            raise unittest.SkipTest("smoke 索引不存在")

    def test_batch_runs_every_record_and_aggregates(self) -> None:
        params = DesignParams(design_scheme="smFISH", desired_probe_count=5)
        seen: list[tuple[int, int, str]] = []
        report = run_batch(
            _two_gene_fasta(),
            str(PROJECT / "data/smoke/target_idx"),
            [],
            params,
            threads=2,
            progress=lambda i, n, rid: seen.append((i, n, rid)),
        )
        self.assertIsInstance(report, BatchReport)
        self.assertEqual(report.total_genes, 2)
        self.assertTrue(all(g.ok for g in report.genes))
        self.assertEqual([g.record_id for g in report.genes], ["geneA", "geneB"])
        # progress 回调按顺序覆盖每条记录
        self.assertEqual([s[0] for s in seen], [1, 2])
        rows = report.per_gene_rows()
        self.assertEqual(len(rows), 2)
        self.assertIn("final_probes", rows[0])
        probe_rows = batch_probe_rows(report)
        for row in probe_rows:
            self.assertIn(row["record_id"], {"geneA", "geneB"})

    def test_single_record_failure_does_not_break_batch(self) -> None:
        bad = ">bad_gene synthetic\nACG\n"  # 太短（<4 nt 会被引擎拒绝）
        good = ">good_gene synthetic\n" + GENE_A + "\n"
        params = DesignParams(design_scheme="smFISH", desired_probe_count=5)
        report = run_batch(
            bad + good, str(PROJECT / "data/smoke/target_idx"), [], params, threads=2
        )
        self.assertFalse(report.genes[0].ok)
        self.assertTrue(report.genes[1].ok)
        self.assertEqual(len(report.failed_genes), 1)
        self.assertGreater(report.total_probes, 0)


if __name__ == "__main__":
    unittest.main()
