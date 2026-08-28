from probedesign.filters import apply_thermo_filters
from probedesign.models import DesignParams, Probe


def make_probe(seq: str) -> Probe:
    return Probe(
        probe_id="test",
        target_id="test",
        start=0,
        stop=len(seq),
        sequence=seq,
        rc_sequence=seq,
    )


def test_gc_filter():
    probe = make_probe("GCGCGCGCGCGCGCGCGCGC")  # 100% GC
    params = DesignParams(min_gc=0.2, max_gc=0.8)
    apply_thermo_filters([probe], params)
    assert not probe.passed
    assert "GC=" in probe.failure_reasons[0]


def test_homopolymer_filter():
    probe = make_probe("ATGCATGCATGCAAAAATGC")
    params = DesignParams(max_homopolymer=4)
    apply_thermo_filters([probe], params)
    assert not probe.passed
    assert any("homopolymer" in r for r in probe.failure_reasons)


def test_passing_probe():
    # Moderate GC, no long homopolymers, reasonable Tm, no stable hairpin.
    # Tm expectation reflects the corrected SantaLucia 1998 NN table.
    probe = make_probe("ATGCGTACGTAGCGAT")
    params = DesignParams(min_tm=45.0, max_tm=70.0)
    apply_thermo_filters([probe], params)
    assert probe.passed
    assert probe.gc_content > 0
    assert probe.tm > 0
