import os
import random
import shutil
import tempfile

import pytest

from mycoprimer.alignment import build_bowtie2_index
from mycoprimer.models import DesignParams, ReferenceGenome
from mycoprimer.pipeline import run_design
from mycoprimer.schemes.hcr3 import design_hcr3
from mycoprimer.schemes.initiators import HCR_INITIATORS
from mycoprimer.schemes.smifish import design_smifish
from mycoprimer.schemes.snail import design_snail


def _random_seq(n: int, gc: float = 0.5, seed: int | None = None) -> str:
    if seed is not None:
        random.seed(seed)
    bases = ["A", "T", "G", "C"]
    weights = [(1 - gc) / 2, (1 - gc) / 2, gc / 2, gc / 2]
    return "".join(random.choices(bases, weights=weights, k=n))


@pytest.fixture
def sample_target():
    tmpdir = tempfile.mkdtemp()
    target_fa = os.path.join(tmpdir, "target.fa")
    host_fa = os.path.join(tmpdir, "host.fa")

    random.seed(42)
    with open(target_fa, "w") as f:
        f.write(">target_1\n" + _random_seq(800, gc=0.5) + "\n")
    with open(host_fa, "w") as f:
        f.write(">host_chr1\n" + _random_seq(3000, gc=0.5) + "\n")

    target_index = os.path.join(tmpdir, "target_idx")
    host_index = os.path.join(tmpdir, "host_idx")
    build_bowtie2_index(target_fa, target_index)
    build_bowtie2_index(host_fa, host_index)

    yield {
        "tmpdir": tmpdir,
        "target_fa": target_fa,
        "host_fa": host_fa,
        "target_index": target_index,
        "host_index": host_index,
    }

    shutil.rmtree(tmpdir, ignore_errors=True)


@pytest.fixture
def host_genome(sample_target):
    return ReferenceGenome(
        id="host",
        organism="Host",
        fasta_path=sample_target["host_fa"],
        bowtie2_index=sample_target["host_index"],
        is_host=True,
    )


def test_smfish_scheme(sample_target, host_genome):
    params = DesignParams(design_scheme="smFISH", desired_probe_count=5)
    result = run_design(
        sample_target["target_fa"],
        sample_target["target_index"],
        [host_genome],
        params,
    )
    assert result.params.design_scheme == "smFISH"
    assert len(result.passed_probes) > 0


def test_smifish_appends_readout(sample_target, host_genome):
    readout = "ATGCATGCATGC"
    params = DesignParams(
        design_scheme="smiFISH",
        desired_probe_count=5,
        smi_readout_sequence=readout,
        smi_readout_position="3prime",
        smi_linker="TTT",
    )
    result = design_smifish(
        sample_target["target_fa"],
        sample_target["target_index"],
        [host_genome],
        params,
    )
    assert len(result.passed_probes) > 0
    for probe in result.passed_probes:
        full = probe.metadata["full_sequence"]
        assert full.endswith(readout)
        assert "TTT" in full
        assert probe.metadata["readout_sequence"] == readout


def test_hcr3_produces_p1_p2(sample_target, host_genome):
    params = DesignParams(
        design_scheme="HCR3",
        hcr_tile_size=52,
        desired_probe_count=5,
        hcr_min_gibbs=-90,
        hcr_max_gibbs=-30,
        hcr_min_gc=30,
        hcr_max_gc=70,
        min_tm=40,
        max_tm=100,
        max_hairpin_tm=70,
        max_homopolymer=7,
        max_target_hits=100,
    )
    result = design_hcr3(
        sample_target["target_fa"],
        sample_target["target_index"],
        [host_genome],
        params,
    )
    assert len(result.passed_probes) > 0
    for probe in result.passed_probes:
        p1 = probe.metadata["P1_sequence"]
        p2 = probe.metadata["P2_sequence"]
        assert HCR_INITIATORS["B1"]["odd"] in p1
        assert HCR_INITIATORS["B1"]["even"] in p2
        assert "gibbs_fe" in probe.metadata
        assert params.hcr_min_gibbs <= probe.metadata["gibbs_fe"] <= params.hcr_max_gibbs


def test_snail_produces_primer_padlock(sample_target, host_genome):
    ugi = "ACGTACGTACGTACGTACGTAC"
    params = DesignParams(
        design_scheme="SNAIL-FISH",
        desired_probe_count=5,
        snail_arm_length=20,
        snail_arm_spacer=1,
        snail_min_gc=30,
        snail_max_gc=70,
        snail_hairpin_dg=-5.0,
        snail_ugi_sequence=ugi,
        max_target_hits=100,
    )
    result = design_snail(
        sample_target["target_fa"],
        sample_target["target_index"],
        [host_genome],
        params,
    )
    assert len(result.passed_probes) > 0
    for probe in result.passed_probes:
        primer = probe.metadata["primer_sequence"]
        padlock = probe.metadata["padlock_sequence"]
        assert len(primer) == params.snail_arm_length + len(params.snail_primer_end)
        assert ugi in padlock
        assert params.snail_padlock_start in padlock
        assert params.snail_padlock_end in padlock
        assert len(probe.metadata["arm1_sequence"]) == params.snail_arm_length
        assert len(probe.metadata["arm2_sequence"]) == params.snail_arm_length
