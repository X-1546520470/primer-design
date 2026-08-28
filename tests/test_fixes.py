"""Regression tests for the engine fixes (v2.0 refactor)."""

from __future__ import annotations

import math
import unittest

from probedesign.alignment import parse_sam_hit_counts
from probedesign.models import DesignParams, Probe
from probedesign.selection import select_non_overlapping
from probedesign.utils import calc_tm


class TmRegressionTests(unittest.TestCase):
    """The old NN table mixed up dinucleotides and skipped AA/TT/GG/CC."""

    def test_tm_matches_primer3_within_tolerance(self) -> None:
        try:
            from primer3 import calc_tm as p3tm
        except ImportError:  # pragma: no cover
            self.skipTest("primer3-py not installed")
        sequences = [
            "GCGCGCTTTTTGCGCGC",
            "ATGACCATGATTACGCCAAG",
            "AAAAAAAAAAAAAAAA",
            "ACGCGT",
            "GATTACAGATTACAGATTAC",
        ]
        for seq in sequences:
            ours = calc_tm(seq)
            ref = p3tm(
                seq,
                mv_conc=390,
                dv_conc=0,
                dntp_conc=0,
                dna_conc=250,
                tm_method="santalucia",
                salt_corrections_method="santalucia",
            )
            self.assertLess(abs(ours - ref), 3.5, seq)

    def test_homopolymer_dinucleotides_use_correct_steps(self) -> None:
        # GC-rich stacks (-10.6/-27.2, -9.8/-24.4, -8.0/-19.9) melt far
        # higher than AT-rich stacks (-7.2/-20.4): sanity-check the table
        # produces the expected GC/AT ordering.
        self.assertGreater(calc_tm("GCGCGCGC"), calc_tm("ATATATAT"))


class SamHitCountTests(unittest.TestCase):
    """The old parser skipped secondary (flag 0x100) records, capping every
    probe at 1 hit and silently disabling repeat/host filtering."""

    def test_counts_secondary_alignments(self) -> None:
        import tempfile
        from pathlib import Path

        sam = "\n".join(
            [
                "1\t0\tchr1\t100\t30\t18M\t*\t0\t0\tACGT\tIIII",
                "1\t256\tchr1\t500\t25\t18M\t*\t0\t0\tACGT\tIIII",
                "1\t256\tchr1\t900\t20\t18M\t*\t0\t0\tACGT\tIIII",
                "2\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII",  # unmapped
            ]
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "x.sam"
            path.write_text(sam + "\n")
            counts = parse_sam_hit_counts(str(path), expected=2)
        self.assertEqual(counts["1"], 3)
        self.assertEqual(counts["2"], 0)


class SelectionOverlapTests(unittest.TestCase):
    """The old code compared start distances, letting overlapping probes through."""

    @staticmethod
    def _probe(pid: str, start: int, stop: int, score: float) -> Probe:
        return Probe(
            probe_id=pid,
            target_id="t",
            start=start,
            stop=stop,
            sequence="A" * (stop - start),
            rc_sequence="T" * (stop - start),
            score=score,
        )

    def test_overlapping_probe_is_rejected(self) -> None:
        high = self._probe("high", 10, 34, score=2.0)  # best scoring, 24 nt
        tail = self._probe("tail", 28, 46, score=1.0)  # overlaps 28-34
        selected = select_non_overlapping([high, tail], min_gap=0)
        self.assertEqual([p.probe_id for p in selected], ["high"])

    def test_gap_respected(self) -> None:
        a = self._probe("a", 0, 18, score=2.0)
        b = self._probe("b", 20, 38, score=1.0)
        selected = select_non_overlapping([a, b], min_gap=5)
        self.assertEqual(len(selected), 1)  # 2 nt apart < 5 nt gap

    def test_adjacent_probes_without_gap_both_kept(self) -> None:
        a = self._probe("a", 0, 18, score=2.0)
        b = self._probe("b", 18, 36, score=1.0)
        selected = select_non_overlapping([a, b], min_gap=0)
        self.assertEqual(len(selected), 2)


class SnailArmStrandTests(unittest.TestCase):
    """Primer/padlock arms must be antisense to bind the target RNA."""

    def test_assembled_oligos_hybridize_to_target(self) -> None:
        from probedesign.schemes.snail import (
            _assemble_oligos,
            _filter_arms,
            _mine_snail_candidates,
        )
        from probedesign.utils import reverse_complement

        # Target region chosen so both arms pass the SNAIL filters
        # (55% GC, no repeat motifs, no stable hairpins).
        arm = "GCGTACACGCGTATATACGC"
        spacer_t = "A"
        target = arm + spacer_t + reverse_complement(arm) + "GGGG" * 0
        params = DesignParams(
            design_scheme="SNAIL-FISH",
            snail_arm_length=20,
            snail_arm_spacer=1,
        )
        candidates = _mine_snail_candidates(target, params)
        self.assertTrue(candidates)
        passing = []
        for probe in candidates:
            if _filter_arms(probe, params):
                _assemble_oligos(probe, params)
                passing.append(probe)
        self.assertTrue(passing)
        probe = passing[0]
        arm1_target = target[probe.start : probe.start + 20]
        arm2_target = target[probe.start + 21 : probe.start + 41]
        # The primer's binding arm must be the reverse complement of the
        # target-sense arm1 region (it has to hybridize to the RNA).
        self.assertEqual(
            probe.metadata["arm1_sequence"], reverse_complement(arm1_target)
        )
        self.assertEqual(
            probe.metadata["arm2_sequence"], reverse_complement(arm2_target)
        )
        self.assertTrue(probe.metadata["primer_sequence"].startswith(
            reverse_complement(arm1_target)
        ))


if __name__ == "__main__":
    unittest.main()
