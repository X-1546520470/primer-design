"""高层入口：按参数中的 design_scheme 分发到对应方案模块。

方案注册表见 schemes/__init__.py；新增方案时在 schemes/ 下新建模块并在
注册表中登记即可，pipeline 与 GUI 无需改动。
"""

from __future__ import annotations

from typing import List

from probedesign.models import DesignParams, DesignResult, ReferenceGenome
from probedesign.schemes import design_for_scheme
from probedesign.schemes.smfish import design_smfish


def run_design(
    target_fasta: str,
    target_index: str,
    host_genomes: List[ReferenceGenome],
    params: DesignParams | None = None,
    threads: int = 1,
) -> DesignResult:
    """执行一次探针设计任务。

    参数：
        target_fasta   靶序列 FASTA 文件路径
        target_index   靶标基因组的 bowtie2 索引前缀
        host_genomes   背景/宿主基因组列表（ReferenceGenome）
        params         设计参数；None 时用 smFISH 默认参数
        threads        bowtie2 线程数
    """
    return design_for_scheme(target_fasta, target_index, host_genomes, params, threads=threads)
