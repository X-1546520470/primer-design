"""默认配置与热力学常数。

Tm 计算采用 SantaLucia 1998 统一最近邻模型（unified NN），即探针的熔解温度
由相邻两个核苷酸的"堆积"贡献逐步累加而成。本模块存放：
    - DNA_NN3            16 种二核苷酸的堆积焓/熵（核心参数表）
    - TERMINAL_AT_PENALTY 末端 A/T 的起始修正
    - INITIATION         双链起始焓/熵
    - DEFAULT_*          默认杂交条件（盐浓度、探针浓度、甲酰胺）

注意：v1 的参数表存在错值（混用互补二核苷酸、缺少 AA/TT/GG/CC 四种堆积），
已于 v2 全部修正；修复后的 Tm 与 primer3 独立实现偏差 < 3 °C。
（文献：SantaLucia J. PNAS 1998;95:1460–1465）
"""

from __future__ import annotations

from mycoprimer.models import DesignParams

DEFAULT_PARAMS = DesignParams()

# SantaLucia 1998 统一最近邻参数表，键为探针链 5'->3' 方向的二核苷酸。
# dH 单位 kcal/mol（负值，放能），dS 单位 cal/(mol*K)（负值）。
# 对称的二核苷酸取同一组值（AA/TT 互补对称、GG/CC 对称）；
# 表必须包含全部 16 种组合，否则 calc_tm 会退化到"平均值"近似。
DNA_NN3 = {
    ("A", "A"): (-7.9, -22.2),
    ("T", "T"): (-7.9, -22.2),
    ("A", "T"): (-7.2, -20.4),
    ("T", "A"): (-7.2, -21.3),
    ("A", "G"): (-7.8, -21.0),  # AG/CT
    ("C", "T"): (-7.8, -21.0),
    ("G", "A"): (-8.2, -22.2),  # GA/CT
    ("T", "C"): (-8.2, -22.2),
    ("G", "G"): (-8.0, -19.9),
    ("C", "C"): (-8.0, -19.9),
    ("G", "C"): (-10.6, -27.2),
    ("C", "G"): (-9.8, -24.4),
    ("A", "C"): (-8.4, -22.4),  # AC/GT
    ("G", "T"): (-8.4, -22.4),
    ("C", "A"): (-8.5, -22.7),  # CA/GT
    ("T", "G"): (-8.5, -22.7),
}

# 末端 A/T 起始罚分：每端若为 A/T 碱基对，则额外加上 (+2.2 kcal/mol, +6.9 cal/mol·K)。
# 与 Biopython DNA_NN3 的 init_A/T 约定一致（注意符号为正——v1 曾误写为负）。
TERMINAL_AT_PENALTY = (2.2, 6.9)

# 双链起始焓/熵（SantaLucia 1998），整体计算的初始项。
INITIATION = (0.2, -5.7)

# Tm 计算的默认杂交条件。
# 0.39 M 单价盐 ≈ 2×SSC（0.3 M NaCl + 0.03 M 柠檬酸钠），即 FISH 标准杂交液。
DEFAULT_NA = 0.39  # M
DEFAULT_MG = 0.0  # M
DEFAULT_DNTP = 0.0  # M
DEFAULT_PROBE_CONC = 1e-6  # M
DEFAULT_FORMAMIDE_PCT = 0.0
DEFAULT_FORMAMIDE_FACTOR = 0.65  # degC per percent formamide
