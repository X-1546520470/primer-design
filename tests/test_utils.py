import math

from probedesign.utils import gc_content, has_homopolymer, reverse_complement, calc_tm


def test_reverse_complement():
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("AAAA") == "TTTT"


def test_gc_content():
    assert gc_content("ATGC") == 0.5
    assert gc_content("AAAA") == 0.0
    assert gc_content("GCGC") == 1.0


def test_has_homopolymer():
    assert has_homopolymer("AAAAAT", max_run=4) is True
    assert has_homopolymer("AAAAT", max_run=4) is False
    assert has_homopolymer("AAAAT", max_run=5) is False
    assert has_homopolymer("ATGCATGC", max_run=4) is False


def test_calc_tm_sanity():
    tm = calc_tm("ATGCATGCATGCATGC")
    assert 0 < tm < 100
    # GC-rich should have higher Tm
    tm_gc = calc_tm("GCGCGCGCGCGCGCGC")
    tm_at = calc_tm("ATATATATATATATAT")
    assert tm_gc > tm_at
