"""smiFISH 探针设计 = smFISH + 共享 readout 延伸段（FLAP）。

smiFISH 的探针本身不带荧光：每条探针末端带一段相同的 readout 延伸段，
荧光二级探针与延伸段杂交发声。好处是换荧光颜色只需换二级探针。

实现上复用 design_smfish 的全部流程，仅在最后给每条探针拼接：
    5′ 端方案：readout + linker + 探针
    3′ 端方案：探针 + linker + readout
拼接结果存入 metadata["full_sequence"]，订购时使用完整序列。
"""

from __future__ import annotations

from typing import List

from mycoprimer.models import DesignParams, DesignResult, ReferenceGenome
from mycoprimer.schemes.smfish import design_smfish


def design_smifish(
    target_fasta: str,
    target_index: str,
    host_genomes: List[ReferenceGenome],
    params: DesignParams | None = None,
    threads: int = 1,
) -> DesignResult:
    """先跑 smFISH 设计，再为所有探针拼接共享 readout 延伸段。"""
    params = params or DesignParams(design_scheme="smiFISH")
    result = design_smfish(target_fasta, target_index, host_genomes, params, threads=threads)

    readout = (params.smi_readout_sequence or "").upper()
    linker = params.smi_linker.upper()
    position = params.smi_readout_position

    # V2 新增：readout/linker 只允许 ACGT，避免把简并码或错误序列拼进
    # 订购 oligo（readout 决定二级探针结合，拼错整批探针报废）。
    for label, seq in (("readout", readout), ("linker", linker)):
        invalid = set(seq) - set("ACGT")
        if invalid:
            raise ValueError(
                f"smiFISH {label} 序列含非法碱基 {''.join(sorted(invalid))}，"
                "仅接受 A/C/G/T。"
            )

    for probe in result.probes:
        full_seq = probe.sequence
        if readout:
            if position == "5prime":
                full_seq = readout + linker + probe.sequence
            else:
                full_seq = probe.sequence + linker + readout
        probe.metadata["full_sequence"] = full_seq
        probe.metadata["readout_sequence"] = readout
        probe.metadata["readout_position"] = position
        probe.metadata["linker"] = linker

    return result
