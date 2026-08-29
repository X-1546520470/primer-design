"""设计方案注册与分发。

每个方案是 schemes/ 下的一个独立模块，实现统一签名的 design_* 函数：
    design_xxx(target_fasta, target_index, host_genomes, params, threads) -> DesignResult

注册表 _DESIGN_FUNCS 把 params.design_scheme 名称映射到实现。**新增方案**：
    1. 新建 schemes/my_method.py，实现上述签名的设计函数；
    2. 在本文件 import 并登记到 _DESIGN_FUNCS；
    3. （可选）在 models.DesignParams 中补充方案专属参数，
       在 GUI 的 SCHEME_INFO 中补一句说明。
引擎其余部分与 GUI 均无需改动。
"""

from __future__ import annotations

from typing import Dict, List

from mycoprimer.models import DesignParams, DesignResult, ReferenceGenome
from mycoprimer.schemes.hcr3 import design_hcr3
from mycoprimer.schemes.smfish import design_smfish
from mycoprimer.schemes.smifish import design_smifish
from mycoprimer.schemes.snail import design_snail

# 方案名 -> 设计函数。GUI 侧栏的方案下拉框与这里的键保持一致。
_DESIGN_FUNCS: Dict[str, callable] = {
    "smFISH": design_smfish,
    "smiFISH": design_smifish,
    "HCR3": design_hcr3,
    "SNAIL-FISH": design_snail,
}


def design_for_scheme(
    target_fasta: str,
    target_index: str,
    host_genomes: List[ReferenceGenome],
    params: DesignParams | None = None,
    threads: int = 1,
) -> DesignResult:
    """Dispatch to the appropriate design scheme."""
    params = params or DesignParams()
    scheme = params.design_scheme
    func = _DESIGN_FUNCS.get(scheme)
    if func is None:
        raise ValueError(f"Unknown design scheme: {scheme}")
    return func(target_fasta, target_index, host_genomes, params, threads=threads)
