"""序列工具：反向互补、GC 分数、同聚碱基检测与最近邻 Tm 计算。

Tm 公式（SantaLucia 1998 统一最近邻模型）：

    ΔH = Σ 堆积焓 + 起始焓 + 末端 A/T 修正
    ΔS = Σ 堆积熵 + 起始熵 + 末端 A/T 修正 + 盐修正（0.368·(N−1)·ln[Na⁺]）
    Tm = 1000·ΔH / (ΔS + R·ln(Ct/4)) − 273.15 − 甲酰胺×0.65

其中 Ct/4 为非自互补双链的浓度项；单价盐有效浓度可用 Mg²⁺/dNTP 做von Ahsen
校正（na + 120·√(mg − dntp)）。参数表见 config.DNA_NN3。
"""

from __future__ import annotations

import math

from mycoprimer.config import (
    DEFAULT_DNTP,
    DEFAULT_FORMAMIDE_FACTOR,
    DEFAULT_FORMAMIDE_PCT,
    DEFAULT_MG,
    DEFAULT_NA,
    DEFAULT_PROBE_CONC,
    DNA_NN3,
    INITIATION,
    TERMINAL_AT_PENALTY,
)

COMPLEMENT = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")

_GAS_CONSTANT = 1.9872041  # cal / (mol * K)


def reverse_complement(seq: str) -> str:
    """返回序列的反向互补链（支持 ACGTUN 与小写）。

    探针序列约定：探针与靶 RNA 反义配对，因此"靶标窗口的正链序列"取反向
    互补后才得到真正合成、能杂交上去的探针序列。
    """
    return seq.translate(COMPLEMENT)[::-1]


def gc_content(seq: str) -> float:
    """返回 GC 分数（0-1 的小数；界面展示时乘 100）。空序列返回 0。"""
    if not seq:
        return 0.0
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq)


def has_homopolymer(seq: str, max_run: int = 4) -> bool:
    """检测是否存在长度超过 max_run 的连续相同碱基（如 GGGGG）。

    长同聚串会引物滑移 / 非特异结合，是探针设计的常规过滤项。
    """
    if max_run <= 0:
        return False
    seq = seq.upper()
    count = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            count += 1
            if count > max_run:
                return True
        else:
            count = 1
    return False


def salt_correction(
    na: float = DEFAULT_NA, mg: float = DEFAULT_MG, dntp: float = DEFAULT_DNTP
) -> float:
    """计算有效单价盐浓度（von Ahsen 2001）：Na + 120·√(Mg²⁺ − dNTP)。

    Mg²⁺ 对双链稳定的贡献可折算为高价单价盐；被 dNTP 螯合的部分不参与。
    """
    return na + 120.0 * math.sqrt(max(0.0, mg - dntp))


def calc_tm(
    seq: str,
    probe_conc: float = DEFAULT_PROBE_CONC,
    na: float = DEFAULT_NA,
    mg: float = DEFAULT_MG,
    dntp: float = DEFAULT_DNTP,
    formamide_pct: float = DEFAULT_FORMAMIDE_PCT,
    formamide_factor: float = DEFAULT_FORMAMIDE_FACTOR,
) -> float:
    """计算熔解温度 Tm（°C）——SantaLucia 1998 统一最近邻模型。

    计算步骤：
        1. 从起始项 (0.2, -5.7) 出发；
        2. 每端为 A/T 时加末端罚分 (+2.2, +6.9)；
        3. 累加全部相邻二核苷酸的堆积焓/熵（config.DNA_NN3）；
        4. 熵做盐修正：ΔS += 0.368·(N−1)·ln[单价盐]；
        5. 浓度项取 ln(Ct/4)，即非自互补双链假设；
        6. 甲酰胺按 0.65 °C/% 线性降低 Tm。
    """
    seq = seq.upper()
    if len(seq) == 0:
        return 0.0

    dh, ds = INITIATION

    # terminal AT penalty (positive contribution per SantaLucia/Biopython)
    if seq[0] in ("A", "T"):
        dh += TERMINAL_AT_PENALTY[0]
        ds += TERMINAL_AT_PENALTY[1]
    if seq[-1] in ("A", "T"):
        dh += TERMINAL_AT_PENALTY[0]
        ds += TERMINAL_AT_PENALTY[1]

    for i in range(len(seq) - 1):
        pair = (seq[i], seq[i + 1])
        if pair in DNA_NN3:
            step_dh, step_ds = DNA_NN3[pair]
            dh += step_dh
            ds += step_ds
        else:
            # ambiguous bases: mean stacking values
            dh += -8.0
            ds += -21.0

    monovalent = salt_correction(na, mg, dntp)
    ds += 0.368 * (len(seq) - 1) * math.log(monovalent)

    # Non-self-complementary concentration term: Ct/4.
    tm = (dh * 1000.0) / (ds + _GAS_CONSTANT * math.log(probe_conc / 4.0)) - 273.15
    tm -= formamide_pct * formamide_factor
    return tm


def calc_tm_batch(sequences: dict[str, str], **kwargs) -> dict[str, float]:
    """Calculate Tm for many sequences."""
    return {name: calc_tm(seq, **kwargs) for name, seq in sequences.items()}
