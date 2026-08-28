"""ProbeStudio — FISH probe design GUI (smFISH / smiFISH / HCR 3.0 / SNAIL)."""

from __future__ import annotations

import html as html_module
import json
from dataclasses import asdict
from pathlib import Path

import pandas as pd
import streamlit as st

from probedesign import __version__
from probedesign.alignment import AlignmentError, build_bowtie2_index
from probedesign.models import DesignParams, DesignResult, ReferenceGenome
from probedesign.pipeline import run_design
from probedesign.report import probes_to_dataframe

APP_DIR = Path(__file__).resolve().parent
REGISTRY_PATH = APP_DIR / "genome_registry.json"
GENOME_DIR = APP_DIR / "genomes"
INDICES_DIR = APP_DIR / "indices"

SCHEME_INFO = {
    "smFISH": {
        "label": "smFISH",
        "description": "18–24 nt 反义寡核苷酸阵列，每条探针 5′或3′端带荧光基团。"
        "输出探针序列即订购序列。",
    },
    "smiFISH": {
        "label": "smiFISH",
        "description": "smFISH 探针 + 共享 readout 延伸段（FLAP），二级探针结合延伸段发出荧光。"
        "输出含完整序列（探针 + linker + readout）。",
    },
    "HCR3": {
        "label": "HCR 3.0",
        "description": "~52 nt 靶标 tile 拆分为两条 25-mer 半探针，各带分裂 initiator（B1–B5）"
        "以启动杂交链式反应扩增。输出 P1/P2 半探针对。",
    },
    "SNAIL-FISH": {
        "label": "SNAIL FISH",
        "description": "相邻两条 20 nt 靶结合臂分别装配 primer（臂1 + 3′ linker）和"
        " 5′-磷酸化 padlock（5′ anchor + 臂2 + spacer + UGI 条码 + spacer + 3′ anchor），"
        "连接后在靶标上滚环扩增。输出 primer 与 padlock 两条订购序列。",
    },
}

APP_CSS = """
<style>
  .block-container {padding-top: 1.6rem; padding-bottom: 3rem; max-width: 1500px;}
  h1 {letter-spacing: -0.02em;}
  [data-testid="stMetric"] {
    background: #FFFFFF; border: 1px solid #E3EAE6; border-radius: 14px;
    padding: 14px 18px; box-shadow: 0 1px 2px rgba(23,32,38,.05);
  }
  [data-testid="stMetricLabel"] p {font-size: .86rem; color: #5C6B73; font-weight: 600;}
  [data-baseweb="tab-list"] {gap: 2px; border-bottom: 2px solid #E3EAE6;}
  [data-baseweb="tab"] {padding: 10px 16px; font-weight: 600; color: #5C6B73;}
  [data-baseweb="tab"][aria-selected="true"] {color: #6A3AA0;}
  [data-baseweb="tab-highlight"] {background-color: #6A3AA0; height: 3px;}
  .stButton > button {border-radius: 10px;}
  [data-testid="stSidebar"] {background: linear-gradient(180deg, #F7F5FB 0%, #FCFBFE 100%);}
  .pf-kicker {color: #6A3AA0; font-weight: 700; letter-spacing: .1em;
    text-transform: uppercase; font-size: .78rem; margin-bottom: .1rem;}
  .pf-hero-sub {color: #52606D; font-size: 1.02rem; margin-top: -.5rem;}
  .pf-chip {background: #EFE9F7; color: #5B3A8E; border: 1px solid #D8CCEF;
    border-radius: 999px; padding: 2px 10px; font-size: .8rem; font-weight: 600;}
  .pf-chip-row {display: flex; flex-wrap: wrap; gap: 6px; margin-top: 2px;}
  .pf-section-h {font-size: .78rem; font-weight: 700; color: #6A3AA0;
    letter-spacing: .06em; margin: 1rem 0 .2rem 0;}
  .pf-caption {color: #6B7A85; font-size: .8rem; margin-top: -.4rem;}
  .pf-interpret {background: #F7F5FB; border: 1px solid #E4DEF0;
    border-left: 4px solid #6A3AA0; border-radius: 10px; padding: 12px 16px; margin: 4px 0;}
  .pf-interpret-title {font-weight: 700; color: #5B3A8E; margin-bottom: 4px;}
  .pf-interpret ul {margin: 0 0 0 1.1rem; padding: 0;}
  .pf-interpret li {margin: 3px 0; color: #33414B; font-size: .92rem;}
  .pf-mono {font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
    overflow-wrap: anywhere;}
</style>
"""


st.set_page_config(
    page_title="ProbeStudio · FISH 探针设计",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="auto",
)
st.markdown(APP_CSS, unsafe_allow_html=True)


# ---------------------------------------------------------------------------
# Genome registry persistence
# ---------------------------------------------------------------------------

def load_registry() -> dict[str, dict[str, str]]:
    if REGISTRY_PATH.is_file():
        try:
            return json.loads(REGISTRY_PATH.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            return {}
    return {}


def save_registry(registry: dict[str, dict[str, str]]) -> None:
    REGISTRY_PATH.write_text(
        json.dumps(registry, ensure_ascii=False, indent=2), encoding="utf-8"
    )


# ---------------------------------------------------------------------------
# Parameter widgets per scheme
# ---------------------------------------------------------------------------

def render_common_params(params: dict) -> None:
    c1, c2 = st.columns(2)
    params["min_length"] = c1.number_input(
        "最短探针长度 (nt)", 15, 60, params["min_length"], help="候选窗口的最小长度。"
    )
    params["max_length"] = c2.number_input(
        "最长探针长度 (nt)", 15, 60, params["max_length"], help="候选窗口的最大长度。"
    )
    c1, c2 = st.columns(2)
    params["min_tm"] = c1.number_input(
        "Tm 下限 (°C)", 30.0, 90.0, params["min_tm"], 0.5,
        help="低于该 Tm 的候选被过滤。Tm 由 SantaLucia 1998 最近邻模型计算，含单价盐与甲酰胺校正。",
    )
    params["max_tm"] = c2.number_input(
        "Tm 上限 (°C)", 30.0, 95.0, params["max_tm"], 0.5,
        help="高于该 Tm 的候选被过滤。窗口越窄探针均一性越好，但覆盖率越低。",
    )
    c1, c2 = st.columns(2)
    params["min_gc"] = c1.number_input(
        "GC 下限 (0-1)", 0.0, 1.0, params["min_gc"], 0.05,
        help="GC 分数下限。常见参考 0.40–0.60，FISH 探针范围通常更宽。",
    )
    params["max_gc"] = c2.number_input(
        "GC 上限 (0-1)", 0.0, 1.0, params["max_gc"], 0.05,
        help="GC 分数上限。过高 GC 增加 non-specific 结合风险。",
    )


def render_smifish_params(params: dict) -> None:
    render_common_params(params)
    params["target_tm"] = st.number_input(
        "目标 Tm (°C)", 0.0, 95.0, params["target_tm"] or 60.0, 0.5,
        help="打分时优先选择接近该 Tm 的探针；设为 0 关闭偏好。",
    )
    params["max_hairpin_tm"] = st.number_input(
        "发卡 Tm 上限 (°C)", 25.0, 80.0, params["max_hairpin_tm"], 1.0,
        help="primer3 计算的发卡 Tm 超过该值的候选被过滤。"
    )
    params["max_homopolymer"] = st.number_input(
        "最长同聚碱基", 3, 8, params["max_homopolymer"],
        help="超过该长度的连续相同碱基（如 GGGG）会被过滤。",
    )


def render_smifish_scheme_params(params: dict) -> None:
    c1, c2 = st.columns(2)
    params["smi_readout_position"] = c1.selectbox(
        "readout 位置", ("3prime", "5prime"),
        help="延伸段加在探针的 3′端还是 5′端（FLAP 方向）。",
    )
    params["smi_linker"] = c2.text_input(
        "linker 序列", params["smi_linker"], help="探针与 readout 延伸段之间的间隔序列。"
    )
    params["smi_readout_sequence"] = st.text_input(
        "readout 序列（可留空）", params["smi_readout_sequence"] or "",
        help="共享的 FLAP 序列，与荧光二级探针互补。留空则仅输出靶结合探针。"


    )


def render_hcr_params(params: dict) -> None:
    params["hcr_tile_size"] = st.number_input(
        "tile 长度 (nt)", 40, 64, params["hcr_tile_size"],
        help="HCR 3.0 标准 tile 为 52 nt，拆分为两条 25-mer 半探针（去掉中间 2 nt）。"
    )
    params["hcr_channel"] = st.selectbox(
        "initiator 通道", ("B1", "B2", "B3", "B4", "B5"),
        help="分裂 initiator 选择；多色实验为不同靶标选不同通道。",
    )
    c1, c2 = st.columns(2)
    params["hcr_min_gibbs"] = c1.number_input(
        "Gibbs 自由能上限 (kcal/mol)", -90.0, 0.0, params["hcr_min_gibbs"], 1.0,
        help="RNA/DNA 杂交体 Gibbs 自由能窗口下限（Sugimoto 1995 参数）。"
    )
    params["hcr_max_gibbs"] = c2.number_input(
        "Gibbs 自由能下限 (kcal/mol)", -90.0, 0.0, params["hcr_max_gibbs"], 1.0,
        help="窗口上限；官方推荐约 −50 到 −70 kcal/mol。"
    )
    c1, c2 = st.columns(2)
    params["hcr_min_gc"] = c1.number_input("tile GC 下限 (%)", 0.0, 100.0, params["hcr_min_gc"])
    params["hcr_max_gc"] = c2.number_input("tile GC 上限 (%)", 0.0, 100.0, params["hcr_max_gc"])
    params["hcr_dtm_max"] = st.number_input(
        "半探针 dTm 上限 (°C)", 0.0, 20.0, params["hcr_dtm_max"] or 5.0, 0.5,
        help="两条半探针 Tm 之差的最大允许值；过大使两条半探针结合不同步。"
    )
    st.caption(
        "52-mer tile 的整体 Tm 通常远高于 smFISH 窗口，因此 tile Tm 过滤默认关闭，"
        "协议以半探针 dTm + Gibbs 窗口为准。"
    )
    with st.expander("tile Tm 窗口（可选，默认关闭）"):
        c1, c2 = st.columns(2)
        hcr_min_tm = c1.number_input("tile Tm 下限（0 = 关闭）", 0.0, 120.0, 0.0, 1.0)
        hcr_max_tm = c2.number_input("tile Tm 上限（0 = 关闭）", 0.0, 120.0, 0.0, 1.0)
        params["hcr_min_tm"] = hcr_min_tm if hcr_min_tm > 0 else None
        params["hcr_max_tm"] = hcr_max_tm if hcr_max_tm > 0 else None


def render_snail_params(params: dict) -> None:
    c1, c2 = st.columns(2)
    params["snail_arm_length"] = c1.number_input(
        "臂长 (nt)", 15, 30, params["snail_arm_length"], help="每条靶结合臂的长度，标准 20 nt。"
    )
    params["snail_arm_spacer"] = c2.number_input(
        "臂间隔 (nt)", 0, 10, params["snail_arm_spacer"], help="两条靶结合臂之间在靶标上的间隔。"
    )
    c1, c2 = st.columns(2)
    params["snail_min_gc"] = c1.number_input("臂 GC 下限 (%)", 0.0, 100.0, params["snail_min_gc"])
    params["snail_max_gc"] = c2.number_input("臂 GC 上限 (%)", 0.0, 100.0, params["snail_max_gc"])
    params["snail_hairpin_dg"] = st.number_input(
        "臂发卡 dG 阈值 (kcal/mol)", -20.0, 0.0, params["snail_hairpin_dg"], 0.5,
        help="任一臂的发卡自由能低于该值（更负）即过滤；标准 −9 kcal/mol。"
    )
    params["snail_ugi_sequence"] = st.text_input(
        "UGI 条码序列（留空用 N）", params["snail_ugi_sequence"] or "",
        help="正交条码序列；留空则 padlock 中以 N 占位，可后续再填。",
    )
    with st.expander("接头与 anchor 序列（默认来自 SNAIL 标准结构）"):
        params["snail_primer_end"] = st.text_input("primer 3′ linker", params["snail_primer_end"])
        params["snail_padlock_start"] = st.text_input("padlock 5′ anchor", params["snail_padlock_start"])
        params["snail_padlock_end"] = st.text_input("padlock 3′ anchor", params["snail_padlock_end"])
        c1, c2 = st.columns(2)
        params["snail_spacer1"] = c1.text_input("spacer1（臂2 与 UGI 之间）", params["snail_spacer1"])
        params["snail_spacer2"] = c2.text_input("spacer2（UGI 与 3′ anchor 之间）", params["snail_spacer2"])


def render_specificity_params(params: dict) -> None:
    params["max_target_hits"] = st.number_input(
        "靶标基因组最大比对数", 1, 200, params["max_target_hits"],
        help="超过该比对数的候选视为重复序列并过滤（按全部 reported alignments 计数）。"
    )
    params["max_host_hits"] = st.number_input(
        "宿主基因组最大比对数", 0, 50, params["max_host_hits"],
        help="宿主/背景基因组上出现超过该比对数即淘汰；0 表示任何比对都淘汰。"
    )
    with st.expander("bowtie2 高级参数"):
        params["bowtie2_preset"] = st.selectbox(
            "灵敏度预设",
            ("--very-sensitive-local", "--sensitive-local", "--local"),
            help="越敏感越慢；特异性分析建议 very-sensitive-local。",
        )
        params["bowtie2_score_min"] = st.text_input(
            "--score-min", params["bowtie2_score_min"],
            help="G,20,8 表示匹配 ≥20 分；数值越高对错配越严格。",
        )
        params["threads"] = st.number_input("线程数", 1, 16, 2)


def render_condition_params(params: dict) -> None:
    c1, c2, c3 = st.columns(3)
    params["na"] = c1.number_input(
        "Na⁺ (M)", 0.0, 2.0, params["na"], 0.01,
        help="单价盐浓度；2×SSC ≈ 0.39 M Na⁺ 等效。"
    )
    params["mg"] = c2.number_input("Mg²⁺ (M)", 0.0, 0.5, params["mg"], 0.001)
    params["dntp"] = c3.number_input("dNTP (M)", 0.0, 0.05, params["dntp"], 0.001)
    c1, c2 = st.columns(2)
    params["probe_conc"] = c1.number_input(
        "探针浓度 (nM)", 0.01, 10000.0, params["probe_conc"] * 1e9, 1.0,
    ) / 1e9
    params["formamide_pct"] = c2.number_input(
        "甲酰胺 (% v/v)", 0, 80, params["formamide_pct"],
        help="甲酰胺降低 Tm（默认 0.65 °C/%）。杂交温度策略在此体现。",
    )


# ---------------------------------------------------------------------------
# Design run
# ---------------------------------------------------------------------------

def collect_params(ui: dict) -> DesignParams:
    return DesignParams(
        min_length=int(ui["min_length"]),
        max_length=int(ui["max_length"]),
        min_tm=float(ui["min_tm"]),
        max_tm=float(ui["max_tm"]),
        target_tm=float(ui["target_tm"]) if ui.get("target_tm") else None,
        min_gc=float(ui["min_gc"]),
        max_gc=float(ui["max_gc"]),
        max_homopolymer=int(ui["max_homopolymer"]),
        max_hairpin_tm=float(ui["max_hairpin_tm"]),
        bowtie2_preset=ui["bowtie2_preset"],
        bowtie2_score_min=ui["bowtie2_score_min"],
        max_target_hits=int(ui["max_target_hits"]),
        max_host_hits=int(ui["max_host_hits"]),
        min_gap=int(ui["min_gap"]),
        desired_probe_count=int(ui["desired_probe_count"]) if ui["desired_probe_count"] else None,
        strand=ui["strand"],
        design_scheme=ui["design_scheme"],
        smi_readout_sequence=(ui.get("smi_readout_sequence") or "").strip() or None,
        smi_readout_position=ui.get("smi_readout_position", "3prime"),
        smi_linker=ui.get("smi_linker", "TTT"),
        hcr_tile_size=int(ui["hcr_tile_size"]),
        hcr_channel=ui["hcr_channel"],
        hcr_min_gibbs=float(ui["hcr_min_gibbs"]),
        hcr_max_gibbs=float(ui["hcr_max_gibbs"]),
        hcr_dtm_max=float(ui["hcr_dtm_max"]) if ui["hcr_dtm_max"] else None,
        hcr_min_tm=ui.get("hcr_min_tm"),
        hcr_max_tm=ui.get("hcr_max_tm"),
        hcr_min_gc=float(ui["hcr_min_gc"]),
        hcr_max_gc=float(ui["hcr_max_gc"]),
        snail_arm_length=int(ui["snail_arm_length"]),
        snail_arm_spacer=int(ui["snail_arm_spacer"]),
        snail_min_gc=float(ui["snail_min_gc"]),
        snail_max_gc=float(ui["snail_max_gc"]),
        snail_hairpin_dg=float(ui["snail_hairpin_dg"]),
        snail_primer_end=ui["snail_primer_end"],
        snail_padlock_start=ui["snail_padlock_start"],
        snail_padlock_end=ui["snail_padlock_end"],
        snail_spacer1=ui["snail_spacer1"],
        snail_spacer2=ui["snail_spacer2"],
        snail_ugi_sequence=(ui.get("snail_ugi_sequence") or "").strip() or None,
    )


def write_target_fasta(text: str) -> Path:
    GENOME_DIR.mkdir(parents=True, exist_ok=True)
    path = GENOME_DIR / "target_query.fa"
    path.write_text(text.strip() + "\n", encoding="utf-8")
    return path


def build_target_index(fasta_text: str, index_name: str) -> str:
    INDICES_DIR.mkdir(parents=True, exist_ok=True)
    fasta_path = GENOME_DIR / f"{index_name}.fa"
    GENOME_DIR.mkdir(parents=True, exist_ok=True)
    fasta_path.write_text(fasta_text.strip() + "\n", encoding="utf-8")
    prefix = str(INDICES_DIR / index_name)
    build_bowtie2_index(str(fasta_path), prefix, threads=2)
    return prefix


def run_design_job(
    target_fasta_path: Path, target_index: str, hosts: list[ReferenceGenome], params: DesignParams
) -> DesignResult:
    return run_design(str(target_fasta_path), target_index, hosts, params, threads=2)


# ---------------------------------------------------------------------------
# Result rendering
# ---------------------------------------------------------------------------

RESULT_COLUMN_HELP = {
    "probe_id": "探针标识：靶标名:起始-终止（0-based，终止不计入）。",
    "start": "靶标上结合区起点（0-based）。",
    "stop": "靶标上结合区终点（0-based，不含）。",
    "length": "结合区长度。",
    "sequence": "探针序列（反义，5′→3′），可直接与靶 RNA 杂交。",
    "gc_content": "GC 分数（0-1）。",
    "tm": "SantaLucia 1998 最近邻 Tm（°C），含盐与甲酰胺校正。",
    "hairpin_tm": "primer3 预测的发卡 Tm（°C）。",
    "target_hits": "在靶标基因组索引上的比对数（含 secondary 比对）。",
    "host_hits": "各宿主/背景基因组上的比对数。",
    "score": "排序分：特异性越高、Tm 越接近目标、GC 越均衡分越高（描述性，非合格判定）。",
    "passed": "是否通过全部过滤并被选入最终集合。",
    "failure_reasons": "未通过的原因（逐条列出）。",
}


def render_funnel(result: DesignResult) -> None:
    import plotly.graph_objects as go

    THERMO_KEYS = ("GC=", "Tm=", "homopolymer", "hairpin", "dTm", "Gibbs=")

    def has_reason(probe: Probe, *keys: str) -> bool:
        return any(key in reason for reason in probe.failure_reasons for key in keys)

    total = len(result.probes)
    thermo_pass = sum(1 for p in result.probes if not has_reason(p, *THERMO_KEYS))
    spec_pass = sum(
        1 for p in result.probes
        if not has_reason(p, *THERMO_KEYS) and not has_reason(p, "hits")
    )
    final = len(result.passed_probes)

    fig = go.Figure(
        go.Bar(
            x=["全部候选", "热力学通过", "特异性通过", "最终集合"],
            y=[total, thermo_pass, spec_pass, final],
            marker_color=["#C9CDD4", "#9C89D9", "#7A5FC7", "#6A3AA0"],
            text=[total, thermo_pass, spec_pass, final],
            textposition="outside",
        )
    )
    fig.update_layout(
        height=290,
        margin=dict(l=10, r=10, t=30, b=10),
        title="探针漏斗：每一步淘汰多少候选",
        yaxis_title="探针数",
    )
    st.plotly_chart(fig, width="stretch", config={"displayModeBar": False})


def render_coverage(result: DesignResult) -> None:
    import plotly.graph_objects as go

    passed = result.passed_probes
    failed = result.failed_probes
    fig = go.Figure()
    if failed:
        fig.add_trace(
            go.Scatter(
                x=[p.start for p in failed],
                y=[p.tm for p in failed],
                mode="markers",
                marker=dict(color="#C9CDD4", size=6),
                name="未入选",
                text=[p.probe_id for p in failed],
                hovertemplate="%{text}<br>Tm=%{y:.1f}°C<extra></extra>",
            )
        )
    if passed:
        fig.add_trace(
            go.Scatter(
                x=[p.start for p in passed],
                y=[p.tm for p in passed],
                mode="markers",
                marker=dict(
                    color=[p.tm for p in passed],
                    colorscale="Viridis",
                    size=9,
                    line=dict(width=1, color="#172026"),
                ),
                name="最终探针",
                customdata=[[p.gc_content, p.length] for p in passed],
                text=[p.sequence for p in passed],
                hovertemplate=(
                    "%{x}-%{x}+%{customdata[1]}<br>Tm=%{y:.1f}°C "
                    "GC=%{customdata[0]:.2f}<br>%{text}<extra></extra>"
                ),
            )
        )
    fig.update_layout(
        height=320,
        margin=dict(l=10, r=10, t=40, b=10),
        title=f"靶标覆盖图（{result.target_id}，{result.target_length} nt）",
        xaxis_title="靶标坐标 (nt)",
        yaxis_title="Tm (°C)",
    )
    st.plotly_chart(fig, width="stretch")


def render_distributions(result: DesignResult) -> None:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    passed = result.passed_probes
    failed = [p for p in result.failed_probes if p.tm > 0]
    fig = make_subplots(rows=1, cols=2, subplot_titles=("Tm 分布", "GC 分布"))
    for values, name, color, col in (
        ([p.tm for p in failed], "未通过", "#C9CDD4", 1),
        ([p.tm for p in passed], "通过", "#6A3AA0", 1),
        ([p.gc_content for p in failed], "未通过", "#C9CDD4", 2),
        ([p.gc_content for p in passed], "通过", "#6A3AA0", 2),
    ):
        if values:
            fig.add_trace(
                go.Histogram(x=values, name=name, marker_color=color, showlegend=col == 1),
                row=1,
                col=col,
            )
    fig.update_layout(height=300, margin=dict(l=10, r=10, t=30, b=10), bargap=0.05)
    st.plotly_chart(fig, width="stretch", config={"displayModeBar": False})


def scheme_columns(scheme: str) -> list[str]:
    if scheme == "smiFISH":
        return ["full_sequence", "readout_sequence", "readout_position"]
    if scheme == "HCR3":
        return ["P1_sequence", "P2_sequence", "channel", "gibbs_fe", "dTm"]
    if scheme == "SNAIL-FISH":
        return ["primer_sequence", "padlock_sequence", "arm1_sequence", "arm2_sequence", "ugi_barcode"]
    return []


def render_results(result: DesignResult) -> None:
    scheme = result.params.design_scheme
    passed = result.passed_probes
    covered = sum(p.stop - p.start for p in passed)
    coverage_pct = covered / result.target_length * 100 if result.target_length else 0.0

    st.subheader("设计结果")
    cols = st.columns(5)
    cols[0].metric(
        "候选探针", len(result.probes), help="在长度窗口内枚举的全部候选窗口数。"
    )
    cols[1].metric(
        "通过全部过滤",
        len(passed),
        help="通过热力学、特异性过滤且被间距规则选中的探针数。",
    )
    cols[2].metric(
        "靶标覆盖率",
        f"{min(coverage_pct, 100):.0f}%",
        help="最终探针的结合区总长占靶标长度的比例（重叠会降低该值）。",
    )
    if passed:
        tm_series = pd.Series([p.tm for p in passed])
        tm_sd = tm_series.std() if len(tm_series) > 1 else 0.0
        cols[3].metric(
            "Tm 均值 ± SD",
            f"{tm_series.mean():.1f} ± {tm_sd:.1f}",
            help="最终集合的 Tm 均一性；SD 越小杂交条件越一致。",
        )
        cols[4].metric(
            "GC 均值",
            f"{pd.Series([p.gc_content for p in passed]).mean() * 100:.0f}%",
            help="最终集合的平均 GC 分数。",
        )

    render_funnel(result)

    tabs = st.tabs(["📋 结果表", "📍 覆盖图", "📊 分布图", "🔍 单条详情", "📦 导出"])
    dataframe = probes_to_dataframe(result)
    base_cols = [
        "probe_id", "start", "stop", "length", "sequence",
        "gc_content", "tm", "target_hits", "passed", "failure_reasons",
    ]
    extra_cols = [c for c in scheme_columns(scheme) if c in dataframe.columns]
    visible_cols = base_cols[:7] + extra_cols + base_cols[7:]
    visible_cols = [c for c in visible_cols if c in dataframe.columns]

    with tabs[0]:
        view = st.radio("显示范围", ("仅最终集合", "全部候选"), horizontal=True)
        shown = dataframe if view == "全部候选" else dataframe[dataframe["passed"]]
        st.dataframe(
            shown[visible_cols],
            hide_index=True,
            width="stretch",
            column_config={
                "sequence": st.column_config.TextColumn(width="large", help=RESULT_COLUMN_HELP["sequence"]),
                "probe_id": st.column_config.TextColumn(help=RESULT_COLUMN_HELP["probe_id"]),
                "failure_reasons": st.column_config.TextColumn(width="large", help=RESULT_COLUMN_HELP["failure_reasons"]),
                "tm": st.column_config.NumberColumn(format="%.1f", help=RESULT_COLUMN_HELP["tm"]),
                "gc_content": st.column_config.NumberColumn(format="%.2f", help=RESULT_COLUMN_HELP["gc_content"]),
                "target_hits": st.column_config.NumberColumn(help=RESULT_COLUMN_HELP["target_hits"]),
                "start": st.column_config.NumberColumn(help=RESULT_COLUMN_HELP["start"]),
                "stop": st.column_config.NumberColumn(help=RESULT_COLUMN_HELP["stop"]),
                "length": st.column_config.NumberColumn(help=RESULT_COLUMN_HELP["length"]),
            },
        )
        st.caption("💡 悬停列名旁的 ⓘ 图标查看每列含义。")

    with tabs[1]:
        render_coverage(result)
        st.caption("灰点＝未入选候选（悬停可看原因相关坐标），彩色点＝最终探针，颜色映射 Tm。")
    with tabs[2]:
        render_distributions(result)
        st.caption("通过（紫）与未通过（灰）候选的 Tm / GC 分布叠加，帮助判断过滤窗口是否设置过严。")
    with tabs[3]:
        if passed:
            options = [p.probe_id for p in passed]
            chosen = st.selectbox("选择探针", options)
            probe = next(p for p in passed if p.probe_id == chosen)
            st.markdown(f"**{html_module.escape(probe.probe_id)}** · {probe.length} nt · "
                        f"Tm {probe.tm:.1f} °C · GC {probe.gc_content:.2f} · "
                        f"靶标比对 {probe.target_hits}")
            for key, label in (
                ("sequence", "探针序列（反义，直接订购）"),
                ("full_sequence", "完整序列（含 readout）"),
                ("P1_sequence", "P1 半探针（5′ initiator + 3′ 半探针）"),
                ("P2_sequence", "P2 半探针（5′ 半探针 + 3′ initiator）"),
                ("primer_sequence", "primer 寡核苷酸"),
                ("padlock_sequence", "padlock 寡核苷酸（订购时加 5′ 磷酸化 /5Phos/）"),
            ):
                value = probe.metadata.get(key)
                if value:
                    st.markdown(f'<div class="pf-caption">{label}</div>', unsafe_allow_html=True)
                    st.code(value, language=None)
            if probe.failure_reasons:
                st.warning("；".join(probe.failure_reasons))
        else:
            st.info("没有通过的探针；请尝试放宽 Tm/GC 窗口或降低特异性阈值。")
    with tabs[4]:
        st.download_button(
            "下载全部候选（CSV）",
            dataframe.to_csv(index=False).encode("utf-8-sig"),
            file_name=f"probedesign_{scheme}_all.csv",
            mime="text/csv",
        )
        passed_df = dataframe[dataframe["passed"]]
        st.download_button(
            "下载最终集合（CSV）",
            passed_df.to_csv(index=False).encode("utf-8-sig"),
            file_name=f"probedesign_{scheme}_passed.csv",
            mime="text/csv",
        )
        order_col = "full_sequence" if "full_sequence" in passed_df.columns else "sequence"
        order_df = passed_df[["probe_id", order_col]].rename(
            columns={"probe_id": "Name", order_col: "Sequence"}
        )
        if scheme == "SNAIL-FISH" and "padlock_sequence" in passed_df.columns:
            primer_rows = order_df.copy()
            primer_rows["Name"] = primer_rows["Name"] + "_primer"
            primer_rows["Sequence"] = passed_df["primer_sequence"].values
            padlock_rows = order_df.copy()
            padlock_rows["Name"] = padlock_rows["Name"] + "_padlock_5Phos"
            padlock_rows["Sequence"] = passed_df["padlock_sequence"].values
            order_df = pd.concat([primer_rows, padlock_rows], ignore_index=True)
        st.download_button(
            "下载订购表（名称 + 序列，兼容 IDT Bulk Input）",
            order_df.to_csv(index=False).encode("utf-8-sig"),
            file_name=f"probedesign_{scheme}_order.csv",
            mime="text/csv",
        )
        summary = {
            "scheme": scheme,
            "target_id": result.target_id,
            "target_length": result.target_length,
            "total_candidates": len(result.probes),
            "passed": len(passed),
            "params": asdict(result.params),
        }
        st.download_button(
            "下载参数与汇总（JSON）",
            json.dumps(summary, ensure_ascii=False, indent=2).encode("utf-8"),
            file_name=f"probedesign_{scheme}_summary.json",
            mime="application/json",
        )


# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------

ui_params: dict = {}
with st.sidebar:
    st.header("🔬 ProbeStudio")
    st.caption(f"v{__version__} · 本地 FISH 探针设计，序列不出机。")

    st.markdown('<div class="pf-section-h">① 靶序列</div>', unsafe_allow_html=True)
    source_mode = st.radio("输入方式", ("粘贴 FASTA", "使用已注册索引"), horizontal=True)
    strand = st.selectbox("靶标链", ("+", "-"), help="探针设计所依据的链；mRNA 通常选 +。")
    ui_params["strand"] = strand

    st.markdown('<div class="pf-section-h">② 设计方案</div>', unsafe_allow_html=True)
    scheme = st.selectbox(
        "方案",
        list(SCHEME_INFO),
        format_func=lambda key: SCHEME_INFO[key]["label"],
        help="各方案为独立设计模块，可扩展新方案。",
    )
    ui_params["design_scheme"] = scheme
    st.caption(SCHEME_INFO[scheme]["description"])

    registry = load_registry()
    target_index = ""
    target_fasta_text = ""
    if source_mode == "粘贴 FASTA":
        target_fasta_text = st.text_area(
            "靶序列 FASTA", height=170,
            value=">target\nATGACCATGATTACGCCAAGCGCGCTTTTTGCGCGCGATTACAGATTACAGATTACGGCCACTACGGCGTACACGCGTATATACGC",
            help="标准 FASTA；多条记录时仅使用第一条。",
        )
    else:
        registered_targets = {
            key: value for key, value in registry.items() if not value.get("is_host")
        } or registry
        pick = st.selectbox("选择靶标索引", list(registered_targets) or ["（无）"])
        if pick in registered_targets:
            target_index = registered_targets[pick]["index_prefix"]
            st.session_state["pf_target_pick"] = pick

    st.markdown('<div class="pf-section-h">③ 方案参数</div>', unsafe_allow_html=True)
    scheme_params = DesignParams(design_scheme=scheme)
    defaults = asdict(scheme_params)
    ui_params.update({key: defaults[key] for key in defaults})

    if scheme in ("smFISH", "smiFISH"):
        render_smifish_params(ui_params)
        if scheme == "smiFISH":
            render_smifish_scheme_params(ui_params)
    elif scheme == "HCR3":
        render_hcr_params(ui_params)
    elif scheme == "SNAIL-FISH":
        render_snail_params(ui_params)

    params_min_gap = st.number_input(
        "最小间隔 (nt)", 0, 200, 0,
        help="相邻探针结合区之间的最小间隔；SNAIL 方案会自动加上双臂跨度。",
    )
    ui_params["min_gap"] = params_min_gap
    desired = st.number_input(
        "目标探针数（0 = 不限制）", 0, 500, 0,
        help="候选过多时按均匀分布降采样到该数量。",
    )
    ui_params["desired_probe_count"] = desired

    with st.expander("特异性过滤"):
        spec_params = DesignParams()
        spec_defaults = asdict(spec_params)
        ui_params["max_target_hits"] = spec_defaults["max_target_hits"]
        ui_params["max_host_hits"] = spec_defaults["max_host_hits"]
        ui_params["bowtie2_preset"] = spec_defaults["bowtie2_preset"]
        ui_params["bowtie2_score_min"] = spec_defaults["bowtie2_score_min"]
        render_specificity_params(ui_params)

    with st.expander("杂交条件（Tm 计算用）"):
        cond_defaults = asdict(DesignParams())
        ui_params["na"] = 0.39
        ui_params["mg"] = 0.0
        ui_params["dntp"] = 0.0
        ui_params["probe_conc"] = 1e-6
        ui_params["formamide_pct"] = 0
        render_condition_params(ui_params)

    st.markdown('<div class="pf-section-h">④ 宿主/背景基因组</div>', unsafe_allow_html=True)
    host_ids = st.multiselect(
        "比对过滤用的基因组索引",
        list(registry),
        format_func=lambda key: f"{registry[key]['organism']}（{key}）",
        help="设计时将候选探针比对到这些基因组，超阈值的探针被淘汰。",
    )

    run_clicked = st.button("开始设计", type="primary", width="stretch")

# ---------------------------------------------------------------------------
# Genome management page area
# ---------------------------------------------------------------------------

st.markdown('<div class="pf-kicker">Local FISH probe design</div>', unsafe_allow_html=True)
st.title("ProbeStudio")
st.markdown(
    '<p class="pf-hero-sub">smFISH / smiFISH / HCR 3.0 / SNAIL 探针设计 · '
    "宿主基因组过滤 · 全程本地运行</p>",
    unsafe_allow_html=True,
)
st.markdown(
    '<div class="pf-chip-row">'
    '<span class="pf-chip">🖥️ 本地分析 · 序列不出机</span>'
    '<span class="pf-chip">🧬 四种方案 · 独立模块</span>'
    '<span class="pf-chip">📖 每项指标附解释</span>'
    "</div>",
    unsafe_allow_html=True,
)

with st.expander("🧰 基因组与索引管理", expanded=False):
    registry = load_registry()
    if registry:
        st.dataframe(
            pd.DataFrame(
                [
                    {
                        "ID": key,
                        "名称": value["organism"],
                        "索引前缀": value["index_prefix"],
                        "宿主/背景": "是" if value.get("is_host") else "否",
                    }
                    for key, value in registry.items()
                ]
            ),
            hide_index=True,
            width="stretch",
        )
    else:
        st.caption("尚未注册任何基因组索引。")
    st.divider()
    build_cols = st.columns([2, 2, 1])
    new_id = build_cols[0].text_input("索引 ID", value="", placeholder="例如 ecoli_MG1655")
    new_name = build_cols[1].text_input("显示名称", value="", placeholder="例如 E. coli K-12")
    uploaded_fa = build_cols[2].file_uploader("上传基因组 FASTA", type=("fa", "fasta", "fna", "txt"))
    build_clicked = st.button("构建索引并注册", width="stretch")
    if build_clicked:
        if not new_id or uploaded_fa is None:
            st.error("请填写索引 ID 并上传 FASTA 文件。")
        else:
            try:
                GENOME_DIR.mkdir(parents=True, exist_ok=True)
                INDICES_DIR.mkdir(parents=True, exist_ok=True)
                fa_path = GENOME_DIR / f"{new_id}.fa"
                fa_path.write_bytes(uploaded_fa.getvalue())
                prefix = str(INDICES_DIR / new_id)
                with st.spinner("正在构建 bowtie2 索引…"):
                    build_bowtie2_index(str(fa_path), prefix, threads=2)
                registry[new_id] = {
                    "organism": new_name or new_id,
                    "index_prefix": prefix,
                    "fasta_path": str(fa_path),
                    "is_host": True,
                }
                save_registry(registry)
                st.success(f"索引 {new_id} 构建并注册成功。")
                st.rerun()
            except AlignmentError as exc:
                st.error(str(exc))

# ---------------------------------------------------------------------------
# Run design
# ---------------------------------------------------------------------------

if run_clicked:
    try:
        params = collect_params(ui_params)
        hosts = [
            ReferenceGenome(
                id=key,
                organism=registry[key]["organism"],
                fasta_path=registry[key].get("fasta_path", ""),
                bowtie2_index=registry[key]["index_prefix"],
                is_host=True,
            )
            for key in host_ids
            if key in registry
        ]
        if source_mode == "粘贴 FASTA":
            if not target_fasta_text.strip():
                st.error("请粘贴靶序列 FASTA。")
                st.stop()
            fasta_path = write_target_fasta(target_fasta_text)
            prefix = build_target_index(fasta_path.read_text(encoding="utf-8"), "target_query")
        else:
            if not target_index:
                st.error("请选择已注册的靶标索引。")
                st.stop()
            prefix = target_index
            fasta_path = Path(prefix + ".fa")
            if not fasta_path.is_file():
                fasta_path = Path("target.fa")

        progress = st.progress(0.0, text="正在枚举候选并计算热力学参数…")
        result = run_design_job(fasta_path, prefix, hosts, params)
        progress.progress(1.0, text="设计完成")
        st.session_state["probedesign_result"] = result
        st.session_state["probedesign_hosts"] = [h.id for h in hosts]
    except AlignmentError as exc:
        st.error(f"比对失败：{exc}")
    except ValueError as exc:
        st.error(f"参数或输入错误：{exc}")

result = st.session_state.get("probedesign_result")
if isinstance(result, DesignResult):
    hosts_used = st.session_state.get("probedesign_hosts", [])
    if hosts_used:
        st.caption(f"已对比宿主/背景基因组：{'、'.join(hosts_used)}")
    render_results(result)

st.divider()
st.caption(
    "探针序列均为反义（与靶 RNA 互补）；Tm 使用 SantaLucia 1998 最近邻模型"
    "（含单价盐、Mg²⁺/dNTP 与甲酰胺校正）；特异性计数包含 bowtie2 的全部 reported 比对。"
    "所有指标为描述性参考，不构成探针“合格/不合格”判定。"
)
