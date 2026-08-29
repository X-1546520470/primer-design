from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from mycoprimer.mining import mine_candidates
from mycoprimer.models import DesignParams
from mycoprimer.utils import reverse_complement


def test_mine_candidates_count():
    target = SeqRecord(Seq("ATGC" * 50), id="test")
    params = DesignParams(min_length=18, max_length=20)
    candidates = mine_candidates(target, params)
    expected = (len(target) - 18 + 1) + (len(target) - 19 + 1) + (len(target) - 20 + 1)
    assert len(candidates) == expected


def test_mine_candidates_sequence_orientation():
    target = SeqRecord(Seq("ATGCCCCCCCCCCTAA"), id="test")
    params = DesignParams(min_length=16, max_length=16)
    candidates = mine_candidates(target, params)
    assert len(candidates) == 1
    assert candidates[0].sequence == reverse_complement(str(target.seq))
