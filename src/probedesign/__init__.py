"""ProbeStudio —— 本地 FISH 探针设计与过滤引擎。

包结构（自底向上）：
    config.py    热力学参数表（SantaLucia 1998 最近邻）与默认杂交条件
    utils.py     序列工具：Tm 计算、GC、反向互补、同聚碱基检测
    models.py    数据模型：DesignParams（参数）/ Probe（探针）/ DesignResult
    mining.py    候选枚举：在靶序列上按长度窗口滑动取出所有候选
    filters.py   热力学过滤：GC / 同聚碱基 / Tm / 发卡（廉价过滤优先）
    alignment.py bowtie2 封装：比对到靶标与背景基因组并统计命中数
    scoring.py   特异性过滤与打分排序
    selection.py 间距筛选与均匀降采样
    report.py    输出 DataFrame / CSV
    schemes/     各设计方案的独立实现：
                   smfish.py    smFISH（18–24 nt 反义寡核苷酸阵列）
                   smifish.py   smiFISH（smFISH + 共享 readout 延伸段）
                   hcr3.py      HCR 3.0（52-mer tile 拆分双半探针）
                   snail.py     SNAIL FISH（primer + padlock 双臂结构）
    pipeline.py  run_design()：按 params.design_scheme 分发到对应方案

典型用法：
    from probedesign.models import DesignParams, ReferenceGenome
    from probedesign.pipeline import run_design

    params = DesignParams(design_scheme="smFISH")
    result = run_design("target.fa", "target_idx", [背景基因组...], params)
    result.passed_probes   # 通过全部过滤的最终探针

所有指标均为描述性参考，不构成探针"合格/不合格"判定。
"""

from .models import DesignParams, DesignResult, Probe, ReferenceGenome
from .pipeline import run_design
from .report import probes_to_dataframe, write_outputs

__version__ = "2.0.0"

__all__ = [
    "DesignParams",
    "DesignResult",
    "Probe",
    "ReferenceGenome",
    "probes_to_dataframe",
    "run_design",
    "write_outputs",
]
