#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mycoprimer_gui.py — MycoPrimerV2 分枝杆菌 FISH 探针设计工具 (桌面图形界面)

运行:  python3 mycoprimer_gui.py     (或双击 run_gui.command)
风格:  与 T7 盘 Analysis Tools/codon-optimizer 相同的 Tkinter 桌面布局
依赖:  界面仅用 Python 标准库 (tkinter)；设计引擎依赖 mycoprimer 包
       (primer3-py / biopython) 与 bowtie2（比对时调用）

功能:
  - 四种 FISH 方案独立设计: smFISH / smiFISH / HCR 3.0 / SNAIL FISH
  - 设计目标预设: 低丰度转录本检测(液培) / 物种区分探针 / 自定义
  - 背景基因组过滤: 比对到已注册的基因组索引, 命中超阈值即淘汰
  - 输出: 结果表(每列含义见使用说明)、单条详情、CSV/订购表导出、
    基因组与索引管理(上传 FASTA → 构建 bowtie2 索引 → 注册)

设计跑在后台线程, 界面不卡顿; 所有指标为描述性参考, 不构成
探针"合格/不合格"判定。序列全程只在本地处理。
"""

from __future__ import annotations

import csv
import json
import os
import queue
import sys
import threading
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

PROJECT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT / "src"))

from mycoprimer import __version__  # noqa: E402
from mycoprimer.alignment import AlignmentError, build_bowtie2_index  # noqa: E402
from mycoprimer.models import DesignParams, DesignResult, ReferenceGenome  # noqa: E402
from mycoprimer.pipeline import run_design  # noqa: E402

APP_TITLE = f"MycoPrimerV2 分枝杆菌 FISH 探针设计工具 v{__version__}"

DATA_DIR = Path(os.environ.get("PROBESTUDIO_HOME", PROJECT))
REGISTRY_PATH = DATA_DIR / "genome_registry.json"
GENOME_DIR = DATA_DIR / "genomes"
INDICES_DIR = DATA_DIR / "indices"

SCHEMES = {
    "smFISH": "smFISH（18–24 nt 反义寡核苷酸，单碱基级覆盖）",
    "smiFISH": "smiFISH（smFISH + 共享 readout 延伸段 FLAP）",
    "HCR3": "HCR 3.0（52-mer tile 拆分双半探针 + initiator）",
    "SNAIL-FISH": "SNAIL FISH（primer + 5′-磷酸化 padlock 双臂）",
}
SCHEME_NOTES = {
    "smFISH": "每条探针即订购序列；低丰度转录本建议堆叠 ≥40 条。",
    "smiFISH": "换荧光颜色只需换二级探针；readout 请填实验所用 LNA 二级探针的互补序列。",
    "HCR3": "窗口已按分枝杆菌高 GC 基因组(65–67%)调优；多色实验为不同靶标选不同通道。",
    "SNAIL-FISH": "padlock 订购需加 5′ 磷酸化；导出表已自动生成 /5Phos/ 变体。",
}
PRESETS = {
    "低丰度转录本检测（液培）": {"desired": 48, "gap": 2, "note": "不开跨物种背景过滤（单菌培养无跨种背景）"},
    "物种区分探针": {"desired": 20, "gap": 0, "note": "请在背景基因组中勾选需要区分的物种"},
    "自定义": {"desired": 0, "gap": 0, "note": "所有参数手动设置"},
}

HELP_TEXT = """\
【使用方法】
1. ① 粘贴或载入靶序列 FASTA（基因/转录本，几百到几千 nt 最合适；方向 5'→3'）。
2. ② 选设计方案与目标预设，调整参数（每项都有提示），在"背景基因组"页勾选
   需要过滤的基因组（低丰度液培检测通常不需要）。
3. 点"▶ 开始设计"，在右侧页签查看结果并导出。

【设计方案说明】
• smFISH: 18–24 nt 反义寡核苷酸阵列，每条探针 5'/3' 端带荧光；低丰度转录本
  建议堆叠 ≥40 条（信号 ∝ 单转录本结合的探针数）。
• smiFISH: 探针末端带共享 readout 延伸段(FLAP)，荧光二级探针与延伸段杂交。
  换颜色只需换二级探针。readout 序列请填二级探针的互补序列。
• HCR 3.0: 52-mer tile 拆成两条 25-mer 半探针，各带分裂 initiator 启动
  杂交链式反应放大。多色实验为不同靶标选不同通道(B1–B5)。
• SNAIL FISH: 相邻双臂分别装入 primer 与 5'-磷酸化 padlock，连接后滚环
  扩增。padlock 订购需加 /5Phos/（导出表已自动生成变体）；UGI 条码
  留空则用 N 占位，下单前替换为实际正交条码。

【参数说明】
• Tm 窗口: SantaLucia 1998 最近邻模型（含 2×SSC 等效单价盐与甲酰胺校正）。
  窗口越窄探针均一性越好、候选越少。18–24 nt 配 50–70 °C 通常合适。
• GC 窗口: 常见参考 0.40–0.60；FISH 探针常放宽到 0.20–0.80。
• 发卡 Tm 上限: primer3 计算；超过阈值的候选自身折叠会抢占杂交。
• 目标探针数: 候选过多时按位置均匀降采样；0 = 不限制。
• 最小间隔: 相邻探针结合区的最小间距；SNAIL 自动加上双臂跨度。
• max_target_hits: 在靶标基因组上的最大比对数，超出视为重复序列。
• 背景基因组: 勾选后候选会比对到这些基因组，任何命中即淘汰
  （max_host_hits=0）。⚠ 只在需要区分物种时勾选——单菌液培检测
  不需要，开了会把可用探针误杀（见测试报告）。

【结果解读】
• 候选/热力/特异/最终: 漏斗各级存活数；淘汰原因在"淘汰原因"列逐条给出。
• 覆盖率: 最终探针结合区总长 / 靶序列长度；低丰度检测越高越好。
• Tm 均值±SD: 均一性指标，SD 小说明可用同一杂交条件。
• HCR 的"最终探针"= 半探针对(P1+P2)；SNAIL = primer+padlock 对。

【数据与方法】
Tm: SantaLucia 1998 PNAS；盐修正: von Ahsen 2001 (Mg²⁺/dNTP)；甲酰胺 0.65°C/%。
Gibbs (HCR): Sugimoto 1995 RNA/DNA 杂交体参数。
比对: bowtie2 --very-sensitive-local --score-min C,36,0（短探针完美匹配
必须计入；含 secondary 比对）。
本工具仅供科研设计参考，不做"合格/不合格"判定；序列全程本地处理。

【快捷键】
Cmd+1..5 切换右侧页签; 结果表中点击任意行查看单条详情。
"""

DEFAULT_FASTA = """>target
ATGACCATGATTACGCCAAGCGCGCTTTTTGCGCGCGATTACAGATTACAGATTACGGCCACTACGGCGTACACGCGTATATACGCCATCATCATCATCATCATGGCTCGAGCACCAACCGGAAC"""


class ProbeGUI:
    def __init__(self, root: tk.Tk):
        self.root = root
        root.title(APP_TITLE)
        root.geometry("1280x860")
        root.minsize(1080, 740)

        self.scheme_vars: dict[str, tk.Variable] = {}
        self.bg_vars: dict[str, tk.BooleanVar] = {}
        self.result: DesignResult | None = None
        self.ui_queue: queue.Queue = queue.Queue()
        self._target_fasta_hash: str | None = None

        self._build_ui()

        # 后台线程 → UI 的消息轮询
        self.root.after(120, self._poll_queue)

    # ------------------------------------------------------------------
    # 界面构建
    # ------------------------------------------------------------------
    def _build_ui(self):
        outer = ttk.Frame(self.root, padding=8)
        outer.pack(fill="both", expand=True)

        # ---- ① 输入区 ----
        inp = ttk.LabelFrame(outer, text=" ① 靶序列输入 (FASTA, 基因/转录本, 5'→3') ", padding=6)
        inp.pack(fill="x")
        bar = ttk.Frame(inp)
        bar.pack(fill="x")
        ttk.Button(bar, text="载入 FASTA 文件", command=self.load_fasta).pack(side="left")
        ttk.Button(bar, text="载入示例", command=self.load_example).pack(side="left", padx=6)
        ttk.Button(bar, text="清空", command=lambda: self._set_text(self.input_txt, "")).pack(side="left")
        self.detect_var = tk.StringVar(value="待输入…")
        ttk.Label(bar, textvariable=self.detect_var, foreground="#555").pack(side="right")

        wrap = ttk.Frame(inp)
        wrap.pack(fill="both", expand=True)
        self.input_txt = tk.Text(wrap, height=6, wrap="none", undo=True,
                                 font=("Menlo", 12), bg="#fbfbfd")
        ys = ttk.Scrollbar(wrap, orient="vertical", command=self.input_txt.yview)
        self.input_txt.configure(yscrollcommand=ys.set)
        self.input_txt.pack(side="left", fill="both", expand=True)
        ys.pack(side="right", fill="y")
        self.input_txt.bind("<<Modified>>", self._on_input_change)

        # ---- 中部: 左参数 / 右结果 ----
        paned = ttk.PanedWindow(outer, orient="horizontal")
        paned.pack(fill="both", expand=True, pady=(8, 0))
        left = ttk.Frame(paned)
        paned.add(left, weight=0)
        self._build_left(left)
        right = ttk.Frame(paned)
        paned.add(right, weight=1)
        self._build_right(right)

    # ---- 左侧: 方案与参数 ----
    def _build_left(self, parent):
        frame = ttk.LabelFrame(parent, text=" ② 设计方案与参数 ", padding=8)
        frame.pack(fill="both", expand=True)

        ttk.Label(frame, text="设计方案").grid(row=0, column=0, sticky="w")
        self.scheme_var = tk.StringVar(value="smFISH")
        cb = ttk.Combobox(frame, textvariable=self.scheme_var, width=40,
                          values=list(SCHEMES), state="readonly")
        cb.grid(row=1, column=0, sticky="we", pady=(0, 2))
        cb.bind("<<ComboboxSelected>>", lambda _e: self._rebuild_scheme_params())
        self.scheme_note = tk.StringVar(value=SCHEME_NOTES["smFISH"])
        ttk.Label(frame, textvariable=self.scheme_note, foreground="#555",
                  wraplength=380, justify="left").grid(row=2, column=0, sticky="we")

        ttk.Separator(frame).grid(row=3, column=0, sticky="we", pady=6)
        ttk.Label(frame, text="设计目标预设").grid(row=4, column=0, sticky="w")
        self.preset_var = tk.StringVar(value="低丰度转录本检测（液培）")
        pb = ttk.Combobox(frame, textvariable=self.preset_var, width=40,
                          values=list(PRESETS), state="readonly")
        pb.grid(row=5, column=0, sticky="we", pady=(0, 2))
        pb.bind("<<ComboboxSelected>>", lambda _e: self._apply_preset())
        self.preset_note = tk.StringVar(value=PRESETS[self.preset_var.get()]["note"])
        ttk.Label(frame, textvariable=self.preset_note, foreground="#0a7d32",
                  wraplength=380, justify="left").grid(row=6, column=0, sticky="we")

        ttk.Separator(frame).grid(row=7, column=0, sticky="we", pady=6)
        # 方案专属参数容器（切换方案时重建）
        self.scheme_frame = ttk.Frame(frame)
        self.scheme_frame.grid(row=8, column=0, sticky="we")
        self._rebuild_scheme_params()

        ttk.Separator(frame).grid(row=9, column=0, sticky="we", pady=6)
        bgf = ttk.LabelFrame(frame, text=" ③ 背景基因组过滤 (可选, 物种区分时使用) ", padding=4)
        bgf.grid(row=10, column=0, sticky="we", pady=(0, 4))
        self.bg_list_frame = ttk.Frame(bgf)
        self.bg_list_frame.pack(fill="x")
        bgbar = ttk.Frame(bgf)
        bgbar.pack(fill="x", pady=(2, 0))
        ttk.Button(bgbar, text="刷新列表", width=10, command=self.refresh_backgrounds).pack(side="left")
        ttk.Button(bgbar, text="在管理页注册…", width=14,
                   command=lambda: self.nb.select(self.genome_tab)).pack(side="left", padx=6)
        self.bg_hint = tk.StringVar(value="")
        ttk.Label(bgf, textvariable=self.bg_hint, foreground="#555",
                  wraplength=360, justify="left").pack(fill="x")

        # ---- 通用参数 + 运行按钮 ----
        ttk.Separator(frame).grid(row=11, column=0, sticky="we", pady=6)
        gen = ttk.Frame(frame)
        gen.grid(row=12, column=0, sticky="we")
        ttk.Label(gen, text="目标探针数(0=不限)").grid(row=0, column=0, sticky="w")
        self.desired_var = tk.IntVar(value=48)
        ttk.Spinbox(gen, from_=0, to=500, width=6, textvariable=self.desired_var).grid(
            row=0, column=1, sticky="w", padx=(4, 12))
        ttk.Label(gen, text="最小间隔(nt)").grid(row=0, column=2, sticky="w")
        self.gap_var = tk.IntVar(value=2)
        ttk.Spinbox(gen, from_=0, to=200, width=5, textvariable=self.gap_var).grid(
            row=0, column=3, sticky="w", padx=4)

        go = ttk.Frame(frame)
        go.grid(row=13, column=0, sticky="we", pady=(8, 0))
        self.run_btn = ttk.Button(go, text="▶ 开始设计", command=self.start_design)
        self.run_btn.pack(fill="x", ipady=4)
        self.status_var = tk.StringVar(value="就绪")
        ttk.Label(frame, textvariable=self.status_var, foreground="#0a7d32",
                  wraplength=380, justify="left").grid(row=14, column=0, sticky="we", pady=(6, 0))

        frame.columnconfigure(0, weight=1)
        self._apply_preset()
        self.refresh_backgrounds()

    def _apply_preset(self):
        """把所选设计目标预设套用到探针数/间隔等参数（仍可手动修改）。"""
        preset = self.preset_var.get()
        cfg = PRESETS.get(preset)
        if cfg is None:
            return
        self.desired_var.set(cfg["desired"])
        self.gap_var.set(cfg["gap"])
        self.preset_note.set(cfg["note"])

    def _rebuild_scheme_params(self):
        """按所选方案重建参数区（每项带提示 tooltip 式说明文字）。"""
        for child in self.scheme_frame.winfo_children():
            child.destroy()
        self.scheme_vars.clear()
        scheme = self.scheme_var.get()
        self.scheme_note.set(SCHEME_NOTES[scheme])
        v: dict[str, tk.Variable] = self.scheme_vars

        def entry(parent, label, var, row, width=10, hint=""):
            ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w")
            e = ttk.Entry(parent, textvariable=var, width=width)
            e.grid(row=row, column=1, sticky="w", padx=(4, 10))
            if hint:
                ttk.Label(parent, text=hint, foreground="#888").grid(row=row, column=2, sticky="w")

        if scheme == "smFISH":
            box = ttk.LabelFrame(self.scheme_frame, text=" smFISH 参数 ", padding=4)
            box.pack(fill="x")
            v["min_length"] = tk.IntVar(value=18); v["max_length"] = tk.IntVar(value=24)
            self._pair_length_inputs(box, 0)
            v["min_tm"] = tk.DoubleVar(value=50.0); v["max_tm"] = tk.DoubleVar(value=70.0)
            self._pair_tm_inputs(box, 1)
            v["min_gc"] = tk.DoubleVar(value=0.20); v["max_gc"] = tk.DoubleVar(value=0.80)
            self._pair_gc_inputs(box, 2)
            v["max_hairpin_tm"] = tk.DoubleVar(value=45.0)
            entry(box, "发卡 Tm 上限 (°C)", v["max_hairpin_tm"], 3)
            v["max_homopolymer"] = tk.IntVar(value=4)
            entry(box, "最长同聚碱基", v["max_homopolymer"], 4)
        elif scheme == "smiFISH":
            box = ttk.LabelFrame(self.scheme_frame, text=" smiFISH 参数 ", padding=4)
            box.pack(fill="x")
            v["smi_readout_sequence"] = tk.StringVar(value="ACGTCGACTATCGAT")
            entry(box, "readout 序列", v["smi_readout_sequence"], 0, 20,
                  "二级探针互补序列, 订购前替换")
            v["smi_linker"] = tk.StringVar(value="TTT")
            entry(box, "linker", v["smi_linker"], 1, 8)
            v["smi_readout_position"] = tk.StringVar(value="3prime")
            ttk.Label(box, text="readout 位置").grid(row=2, column=0, sticky="w")
            ttk.Combobox(box, textvariable=v["smi_readout_position"], width=8,
                         values=("3prime", "5prime"), state="readonly").grid(
                row=2, column=1, sticky="w", padx=4)
        elif scheme == "HCR3":
            box = ttk.LabelFrame(self.scheme_frame, text=" HCR 3.0 参数 (已按分枝杆菌高 GC 调优) ", padding=4)
            box.pack(fill="x")
            v["hcr_channel"] = tk.StringVar(value="B1")
            ttk.Label(box, text="initiator 通道").grid(row=0, column=0, sticky="w")
            ttk.Combobox(box, textvariable=v["hcr_channel"], width=6,
                         values=("B1", "B2", "B3", "B4", "B5"),
                         state="readonly").grid(row=0, column=1, sticky="w", padx=4)
            v["hcr_tile_size"] = tk.IntVar(value=52)
            entry(box, "tile 长度", v["hcr_tile_size"], 1)
            v["hcr_min_gc"] = tk.DoubleVar(value=40.0); v["hcr_max_gc"] = tk.DoubleVar(value=65.0)
            self._pair_gc_inputs(box, row=2, lo_key="hcr_min_gc", hi_key="hcr_max_gc",
                                 label="tile GC 窗口 (%)")
            v["hcr_min_gibbs"] = tk.DoubleVar(value=-75.0); v["hcr_max_gibbs"] = tk.DoubleVar(value=-45.0)
            ttk.Label(box, text="Gibbs 窗口 (kcal/mol)").grid(row=3, column=0, sticky="w")
            f = ttk.Frame(box); f.grid(row=3, column=1, sticky="w")
            ttk.Entry(f, textvariable=v["hcr_min_gibbs"], width=6).pack(side="left")
            ttk.Label(f, text=" ~ ").pack(side="left")
            ttk.Entry(f, textvariable=v["hcr_max_gibbs"], width=6).pack(side="left")
            v["hcr_dtm_max"] = tk.DoubleVar(value=8.0)
            entry(box, "半探针 dTm 上限 (°C)", v["hcr_dtm_max"], 4)
        elif scheme == "SNAIL-FISH":
            box = ttk.LabelFrame(self.scheme_frame, text=" SNAIL FISH 参数 ", padding=4)
            box.pack(fill="x")
            v["snail_arm_length"] = tk.IntVar(value=20)
            entry(box, "臂长 (nt)", v["snail_arm_length"], 0)
            v["snail_arm_spacer"] = tk.IntVar(value=1)
            entry(box, "臂间隔 (nt)", v["snail_arm_spacer"], 1)
            v["snail_min_gc"] = tk.DoubleVar(value=40.0); v["snail_max_gc"] = tk.DoubleVar(value=63.0)
            self._pair_gc_inputs(box, row=2, lo_key="snail_min_gc", hi_key="snail_max_gc",
                                 label="单臂 GC 窗口 (%)")
            v["snail_hairpin_dg"] = tk.DoubleVar(value=-9.0)
            entry(box, "臂发卡 dG 阈值 (kcal/mol)", v["snail_hairpin_dg"], 3,
                  hint="≤ 该值(更负)即淘汰")
            v["snail_ugi_sequence"] = tk.StringVar(value="")
            entry(box, "UGI 条码 (留空=N)", v["snail_ugi_sequence"], 4, 24,
                  "下单前替换为实际条码")

        box = ttk.LabelFrame(self.scheme_frame, text=" 热力学与过滤 (通用) ", padding=4)
        box.pack(fill="x", pady=(6, 0))
        v["max_target_hits"] = tk.IntVar(value=10)
        entry(box, "靶标基因组最大比对数", v["max_target_hits"], 0,
              hint="超出视为重复序列")

    # 成对输入的便捷控件（长度/Tm/GC 等区间参数）
    def _pair_length_inputs(self, parent, row):
        ttk.Label(parent, text="探针长度范围 (nt)").grid(row=row, column=0, sticky="w")
        f = ttk.Frame(parent); f.grid(row=row, column=1, sticky="w")
        ttk.Entry(f, textvariable=self.scheme_vars["min_length"], width=5).pack(side="left")
        ttk.Label(f, text=" – ").pack(side="left")
        ttk.Entry(f, textvariable=self.scheme_vars["max_length"], width=5).pack(side="left")

    def _pair_tm_inputs(self, parent, row=1):
        ttk.Label(parent, text="Tm 窗口 (°C)").grid(row=row, column=0, sticky="w")
        f = ttk.Frame(parent); f.grid(row=row, column=1, sticky="w")
        ttk.Entry(f, textvariable=self.scheme_vars["min_tm"], width=6).pack(side="left")
        ttk.Label(f, text=" – ").pack(side="left")
        ttk.Entry(f, textvariable=self.scheme_vars["max_tm"], width=6).pack(side="left")

    def _pair_gc_inputs(self, parent, row, lo_key="min_gc", hi_key="max_gc", label="GC 窗口 (0-1)"):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w")
        f = ttk.Frame(parent); f.grid(row=row, column=1, sticky="w")
        ttk.Entry(f, textvariable=self.scheme_vars[lo_key], width=6).pack(side="left")
        ttk.Label(f, text=" – ").pack(side="left")
        ttk.Entry(f, textvariable=self.scheme_vars[hi_key], width=6).pack(side="left")

    # ---- 右侧: 结果 Notebook ----
    def _build_right(self, parent):
        self.nb = ttk.Notebook(parent)
        self.nb.pack(fill="both", expand=True)

        # 结果表
        f1 = ttk.Frame(self.nb, padding=4)
        self.nb.add(f1, text=" 结果表 ")
        b1 = ttk.Frame(f1)
        b1.pack(fill="x")
        ttk.Button(b1, text="复制选中探针序列", command=self.copy_selected).pack(side="left")
        ttk.Label(b1, text="点击行查看单条详情（含方案专属寡核苷酸）",
                  foreground="#555").pack(side="left", padx=8)
        self.count_var = tk.StringVar(value="(尚未运行)")
        ttk.Label(b1, textvariable=self.count_var, foreground="#555").pack(side="right")
        wrap = ttk.Frame(f1)
        wrap.pack(fill="both", expand=True)
        cols = ("start", "stop", "length", "sequence", "tm", "gc", "hits", "reason")
        self.tree = ttk.Treeview(wrap, columns=cols, show="headings", selectmode="browse")
        headers = {"start": "起点", "stop": "终点", "length": "长度",
                   "sequence": "序列 (5'→3')", "tm": "Tm", "gc": "GC",
                   "hits": "靶标比对", "reason": "状态"}
        widths = {"start": 70, "stop": 70, "length": 50, "sequence": 260,
                  "tm": 60, "gc": 55, "hits": 70, "reason": 260}
        for c in cols:
            self.tree.heading(c, text=headers[c])
            self.tree.column(c, width=widths[c], anchor="w")
        ysb = ttk.Scrollbar(wrap, orient="vertical", command=self.tree.yview)
        xsb = ttk.Scrollbar(wrap, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=ysb.set, xscrollcommand=xsb.set)
        self.tree.grid(row=0, column=0, sticky="nsew")
        ysb.grid(row=0, column=1, sticky="ns")
        xsb.grid(row=1, column=0, sticky="we")
        wrap.rowconfigure(0, weight=1)
        wrap.columnconfigure(0, weight=1)
        self.tree.bind("<<TreeviewSelect>>", self._on_row_select)

        # 单条详情
        f2 = ttk.Frame(self.nb, padding=4)
        self.nb.add(f2, text=" 单条详情 ")
        wrap2 = ttk.Frame(f2)
        wrap2.pack(fill="both", expand=True)
        self.detail_txt = self._make_text(wrap2)

        # 设计报告
        f3 = ttk.Frame(self.nb, padding=4)
        self.nb.add(f3, text=" 设计报告 ")
        b3 = ttk.Frame(f3)
        b3.pack(fill="x")
        ttk.Button(b3, text="保存设计报告 (txt)", command=self.save_report).pack(side="left")
        wrap3 = ttk.Frame(f3)
        wrap3.pack(fill="both", expand=True)
        self.report_txt = self._make_text(wrap3)

        # 导出
        f4 = ttk.Frame(self.nb, padding=4)
        self.nb.add(f4, text=" 导出 ")
        ttk.Label(f4, text="设计完成后可导出：", padding=6).pack(anchor="w")
        grid = ttk.Frame(f4, padding=8)
        grid.pack(fill="both")
        ttk.Button(grid, text="最终探针 CSV", command=lambda: self.save_csv(True),
                   width=28).grid(row=0, column=0, sticky="we", pady=4)
        ttk.Button(grid, text="全部候选 CSV", command=lambda: self.save_csv(False),
                   width=28).grid(row=1, column=0, sticky="we", pady=4)
        ttk.Button(grid, text="订购表 (IDT 兼容)", command=self.save_order,
                   width=28).grid(row=2, column=0, sticky="we", pady=4)
        ttk.Label(f4, text="SNAIL 订购表自动包含 /5Phos/ 变体；smiFISH 使用完整序列\n"
                           "(探针 + linker + readout)；HCR 3.0 每对拆为 P1/P2 两行。",
                  foreground="#555", padding=6).pack(anchor="w")

        # 基因组与索引
        self.genome_tab = ttk.Frame(self.nb, padding=4)
        self.nb.add(self.genome_tab, text=" 基因组与索引 ")
        gb = ttk.Frame(self.genome_tab)
        gb.pack(fill="x")
        ttk.Button(gb, text="注册基因组 FASTA…", command=self.register_genome).pack(side="left")
        ttk.Button(gb, text="删除选中", command=self.delete_genome).pack(side="left", padx=6)
        ttk.Button(gb, text="刷新", command=self.refresh_backgrounds).pack(side="left")
        self.genome_list = tk.Listbox(self.genome_tab, height=8, font=("Menlo", 12))
        self.genome_list.pack(fill="both", expand=True, pady=(6, 0))
        ttk.Label(self.genome_tab, text=(
            "注册流程：选择目标/背景物种的基因组 FASTA → 自动构建 bowtie2 索引并登记。\n"
            "已在 T7/项目内预注册：MTB H37Rv、BCG Pasteur、M. smegmatis mc² 155。\n"
            "背景基因组只在\"物种区分\"场景勾选；低丰度液培检测无需勾选。"),
            foreground="#555", justify="left").pack(anchor="w", pady=(6, 0))

        # 使用说明
        f5 = ttk.Frame(self.nb, padding=4)
        self.nb.add(f5, text=" 使用说明 ")
        wrap5 = ttk.Frame(f5)
        wrap5.pack(fill="both", expand=True)
        self.help_txt = self._make_text(wrap5)
        self._set_text(self.help_txt, HELP_TEXT)

    def _make_text(self, parent) -> tk.Text:
        txt = tk.Text(parent, wrap="none", font=("Menlo", 12),
                      bg="#ffffff", relief="flat", state="disabled")
        ys = ttk.Scrollbar(parent, orient="vertical", command=txt.yview)
        xs = ttk.Scrollbar(parent, orient="horizontal", command=txt.xview)
        txt.configure(yscrollcommand=ys.set, xscrollcommand=xs.set)
        txt.grid(row=0, column=0, sticky="nsew")
        ys.grid(row=0, column=1, sticky="ns")
        xs.grid(row=1, column=0, sticky="we")
        parent.rowconfigure(0, weight=1)
        parent.columnconfigure(0, weight=1)
        return txt

    def _set_text(self, widget: tk.Text, content: str):
        widget.configure(state="normal")
        widget.delete("1.0", "end")
        widget.insert("1.0", content)
        widget.configure(state="disabled")

    # ------------------------------------------------------------------
    # 输入区行为
    # ------------------------------------------------------------------
    def _on_input_change(self, _event=None):
        text = self.input_txt.get("1.0", "end")
        records = text.count(">")
        seq = "".join(l.strip() for l in text.splitlines() if not l.startswith(">"))
        self.detect_var.set(
            f"检测: {max(records, 1 if seq else 0)} 条记录 · {len(seq)} nt"
            if seq else "待输入…"
        )

    def load_fasta(self):
        path = filedialog.askopenfilename(
            title="选择 FASTA 文件",
            filetypes=(("FASTA", "*.fa *.fasta *.fna *.txt"), ("所有文件", "*.*")),
        )
        if not path:
            return
        text = Path(path).read_text(encoding="utf-8", errors="replace")
        self._set_text(self.input_txt, text)
        self._on_input_change()

    def load_example(self):
        example = next(
            (
                Path(PROJECT / "test_data/lowabundance_results/gene_fastas/msm_MSMEG_0786.fasta")
                .read_text()
                for p in [1]
                if (PROJECT / "test_data/lowabundance_results/gene_fastas/msm_MSMEG_0786.fasta").is_file()
            ),
            DEFAULT_FASTA,
        )
        self._set_text(self.input_txt, example)
        self._on_input_change()

    # ------------------------------------------------------------------
    # 设计运行（后台线程）
    # ------------------------------------------------------------------
    def _collect_params(self) -> DesignParams:
        scheme = self.scheme_var.get()
        v = {k: var.get() for k, var in self.scheme_vars.items()}
        length = None
        if scheme == "smFISH":
            # 长度输入是独立变量
            length = (v.get("min_length"), v.get("max_length"))
        kwargs = dict(
            design_scheme=scheme,
            desired_probe_count=int(self.desired_var.get()) or None,
            min_gap=int(self.gap_var.get()),
            max_target_hits=int(v.get("max_target_hits", 10)),
            strand="+",
        )
        if scheme in ("smFISH", "smiFISH"):
            kwargs.update(
                min_length=int(length[0]) if length else 18,
                max_length=int(length[1]) if length else 24,
                min_tm=float(v.get("min_tm", 50.0)),
                max_tm=float(v.get("max_tm", 70.0)),
                min_gc=float(v.get("min_gc", 0.20)),
                max_gc=float(v.get("max_gc", 0.80)),
                max_hairpin_tm=float(v.get("max_hairpin_tm", 45.0)),
                max_homopolymer=int(v.get("max_homopolymer", 4)),
                target_tm=None,
            )
        if scheme == "smiFISH":
            kwargs.update(
                smi_readout_sequence=v.get("smi_readout_sequence", "").strip() or None,
                smi_readout_position=v.get("smi_readout_position", "3prime"),
                smi_linker=v.get("smi_linker", "TTT"),
            )
        if scheme == "HCR3":
            kwargs.update(
                hcr_channel=v.get("hcr_channel", "B1"),
                hcr_tile_size=int(v.get("hcr_tile_size", 52)),
                hcr_min_gc=float(v.get("hcr_min_gc", 40.0)),
                hcr_max_gc=float(v.get("hcr_max_gc", 65.0)),
                hcr_min_gibbs=float(v.get("hcr_min_gibbs", -75.0)),
                hcr_max_gibbs=float(v.get("hcr_max_gibbs", -45.0)),
                hcr_dtm_max=float(v.get("hcr_dtm_max", 8.0)),
            )
        if scheme == "SNAIL-FISH":
            kwargs.update(
                snail_arm_length=int(v.get("snail_arm_length", 20)),
                snail_arm_spacer=int(v.get("snail_arm_spacer", 1)),
                snail_min_gc=float(v.get("snail_min_gc", 40.0)),
                snail_max_gc=float(v.get("snail_max_gc", 63.0)),
                snail_hairpin_dg=float(v.get("snail_hairpin_dg", -9.0)),
                snail_ugi_sequence=v.get("snail_ugi_sequence", "").strip() or None,
            )
        return DesignParams(**kwargs)

    def start_design(self):
        text = self.input_txt.get("1.0", "end").strip()
        if not text:
            messagebox.showwarning("提示", "请先在 ① 输入靶序列 FASTA。")
            return
        try:
            params = self._collect_params()
        except (TypeError, ValueError) as exc:
            messagebox.showerror("参数错误", f"参数填写有误：{exc}")
            return

        # 背景基因组
        registry = load_registry()
        hosts = [
            ReferenceGenome(
                id=gid,
                organism=registry[gid].get("organism", gid),
                fasta_path=registry[gid].get("fasta_path", ""),
                bowtie2_index=registry[gid].get("index_prefix", ""),
                is_host=True,
            )
            for gid, var in self.bg_vars.items() if var.get() and gid in registry
        ]

        self.run_btn.configure(state="disabled")
        self.status_var.set("正在准备…")
        self.ui_queue.put(("status", "正在写入靶序列并构建索引…"))
        threading.Thread(target=self._design_worker, args=(text, params, hosts),
                         daemon=True).start()

    def _design_worker(self, fasta_text: str, params: DesignParams, hosts):
        q = self.ui_queue
        try:
            GENOME_DIR.mkdir(parents=True, exist_ok=True)
            INDICES_DIR.mkdir(parents=True, exist_ok=True)
            fasta_path = GENOME_DIR / "_target_query.fa"
            fasta_path.write_text(fasta_text.strip() + "\n", encoding="utf-8")

            import hashlib
            digest = hashlib.md5(fasta_text.encode()).hexdigest()
            prefix = str(INDICES_DIR / "_target_query")
            if self._target_fasta_hash != digest:
                q.put(("status", "正在构建靶标 bowtie2 索引…"))
                build_bowtie2_index(str(fasta_path), prefix, threads=2)
                self._target_fasta_hash = digest

            q.put(("status", f"正在设计（{params.design_scheme}，"
                             f"{len(hosts)} 个背景基因组）…"))
            result = run_design(str(fasta_path), prefix, hosts, params, threads=2)
            q.put(("done", result))
        except AlignmentError as exc:
            q.put(("error", f"比对失败：{exc}"))
        except ValueError as exc:
            q.put(("error", f"参数或输入错误：{exc}"))
        except Exception as exc:  # 兜底，避免线程静默失败
            q.put(("error", f"设计失败：{exc}"))

    def _poll_queue(self):
        try:
            while True:
                kind, payload = self.ui_queue.get_nowait()
                if kind == "status":
                    self.status_var.set(payload)
                elif kind == "done":
                    self.result = payload
                    self.run_btn.configure(state="normal")
                    self.status_var.set(
                        f"设计完成：候选 {len(payload.probes)} 条，"
                        f"最终 {len(payload.passed_probes)} 条。"
                    )
                    self._show_result(payload)
                elif kind == "genome_ok":
                    self.run_btn.configure(state="normal")
                    self.status_var.set(f"基因组 {payload} 注册成功。")
                    self.refresh_backgrounds()
                elif kind == "error":
                    self.run_btn.configure(state="normal")
                    self.status_var.set("失败")
                    messagebox.showerror("设计失败", payload)
        except queue.Empty:
            pass
        self.root.after(120, self._poll_queue)

    # ------------------------------------------------------------------
    # 结果展示与导出
    # ------------------------------------------------------------------
    def _show_result(self, result: DesignResult):
        self.tree.delete(*self.tree.get_children())
        scheme = result.params.design_scheme
        for p in result.passed_probes:
            self.tree.insert("", "end", iid=p.probe_id, values=(
                p.start + 1, p.stop, p.length, p.sequence,
                f"{p.tm:.1f}", f"{p.gc_content:.2f}", p.target_hits, "✓ 最终"),
                tags=("pass",))
        for p in result.failed_probes:
            self.tree.insert("", "end", iid=p.probe_id, values=(
                p.start + 1, p.stop, p.length, p.sequence,
                f"{p.tm:.1f}", f"{p.gc_content:.2f}", p.target_hits,
                "; ".join(p.failure_reasons)[:80]), tags=("fail",))
        self.tree.tag_configure("pass", foreground="#0a7d32")
        self.tree.tag_configure("fail", foreground="#9aa")

        covered = sum(p.stop - p.start for p in result.passed_probes)
        cov = covered / result.target_length * 100 if result.target_length else 0
        self.count_var.set(
            f"候选 {len(result.probes)} | 最终 {len(result.passed_probes)} | 覆盖 {cov:.0f}%"
        )

        # 设计报告页
        lines = [
            f"设计方案: {scheme}",
            f"靶标: {result.target_id} ({result.target_length} nt)",
            f"背景基因组: {', '.join(result.host_genome_ids) or '无'}",
            "",
            f"候选总数: {len(result.probes)}",
            f"最终探针: {len(result.passed_probes)}",
            f"覆盖率: {cov:.1f}%",
        ]
        if result.passed_probes:
            tms = [p.tm for p in result.passed_probes]
            lines.append(f"Tm 均值±SD: {statistics_mean(tms):.1f} ± {statistics_sd(tms):.1f} °C")
            lines.append(f"GC 均值: {statistics_mean([p.gc_content for p in result.passed_probes]) * 100:.0f}%")
        lines.append("")
        lines.append("参数:")
        for key, value in sorted(vars(result.params).items()):
            lines.append(f"  {key} = {value}")
        lines.append("")
        lines.append("淘汰原因统计:")
        reasons: dict[str, int] = {}
        for p in result.failed_probes:
            for reason in p.failure_reasons:
                key = reason.split("=")[0].split(" ")[0]
                reasons[key] = reasons.get(key, 0) + 1
        for key, count in sorted(reasons.items(), key=lambda kv: -kv[1]):
            lines.append(f"  {key}: {count} 条候选")
        self._set_text(self.report_txt, "\n".join(lines))

    def _on_row_select(self, _event=None):
        if self.result is None:
            return
        selection = self.tree.selection()
        if not selection:
            return
        probe_id = selection[0]
        probe = next((p for p in self.result.probes if p.probe_id == probe_id), None)
        if probe is None:
            return
        lines = [
            f"探针 ID: {probe.probe_id}",
            f"位置: {probe.start + 1}–{probe.stop} (1-based) | 长度 {probe.length} nt",
            f"序列 (反义, 订购用): {probe.sequence}",
            f"Tm: {probe.tm:.1f} °C | GC: {probe.gc_content:.2f} | "
            f"发卡 Tm: {probe.hairpin_tm:.1f} °C",
            f"靶标比对数: {probe.target_hits} | 背景比对: {probe.host_hits or '未启用'}",
            f"打分: {probe.score:.3f}",
            "",
        ]
        meta_order = [
            ("full_sequence", "完整序列 (探针+linker+readout)"),
            ("P1_sequence", "P1 半探针 (5′ initiator + 3′ 半探针)"),
            ("P2_sequence", "P2 半探针 (5′ 半探针 + 3′ initiator)"),
            ("primer_sequence", "primer 寡核苷酸"),
            ("padlock_sequence_5phos", "padlock 寡核苷酸 (含 /5Phos/, 直接订购)"),
            ("arm1_sequence", "结合臂 1 (反义)"),
            ("arm2_sequence", "结合臂 2 (反义)"),
        ]
        for key, label in meta_order:
            value = probe.metadata.get(key)
            if value:
                lines.append(f"{label}: {value}")
        if probe.failure_reasons:
            lines.append("")
            lines.append("淘汰原因: " + "; ".join(probe.failure_reasons))
        self._set_text(self.detail_txt, "\n".join(lines))

    def copy_selected(self):
        selection = self.tree.selection()
        if not selection:
            messagebox.showinfo("提示", "请先在结果表中选中一条探针。")
            return
        probe = next((p for p in self.result.probes if p.probe_id == selection[0]), None)
        if probe:
            self.root.clipboard_clear()
            self.root.clipboard_append(probe.sequence)
            self.status_var.set(f"已复制 {probe.probe_id} 的序列。")

    # ---- 导出 ----
    def _probe_rows(self, result: DesignResult, passed_only: bool):
        scheme = result.params.design_scheme
        rows = []
        for p in result.probes:
            if passed_only and not p.passed:
                continue
            row = {
                "probe_id": p.probe_id, "start": p.start + 1, "stop": p.stop,
                "length": p.length, "sequence": p.sequence,
                "tm": round(p.tm, 2), "gc": round(p.gc_content, 3),
                "target_hits": p.target_hits, "passed": p.passed,
                "failure_reasons": "; ".join(p.failure_reasons),
            }
            for key in ("full_sequence", "P1_sequence", "P2_sequence",
                        "primer_sequence", "padlock_sequence_5phos"):
                if key in p.metadata:
                    row[key] = p.metadata[key]
            rows.append(row)
        return rows

    def save_csv(self, passed_only: bool):
        if self.result is None:
            messagebox.showinfo("提示", "请先运行一次设计。")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            initialfile=f"probes_{'passed' if passed_only else 'all'}.csv",
            filetypes=(("CSV", "*.csv"),),
        )
        if not path:
            return
        rows = self._probe_rows(self.result, passed_only)
        with open(path, "w", newline="", encoding="utf-8-sig") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        self.status_var.set(f"已导出 {len(rows)} 行 → {path}")

    def save_order(self):
        if self.result is None:
            messagebox.showinfo("提示", "请先运行一次设计。")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".csv", initialfile="order_sheet.csv",
            filetypes=(("CSV", "*.csv"),),
        )
        if not path:
            return
        scheme = self.result.params.design_scheme
        out = []
        for p in self.result.passed_probes:
            if scheme == "smiFISH" and p.metadata.get("full_sequence"):
                out.append({"Name": p.probe_id, "Sequence": p.metadata["full_sequence"]})
            elif scheme == "HCR3":
                if p.metadata.get("P1_sequence"):
                    out.append({"Name": f"{p.probe_id}_P1", "Sequence": p.metadata["P1_sequence"]})
                if p.metadata.get("P2_sequence"):
                    out.append({"Name": f"{p.probe_id}_P2", "Sequence": p.metadata["P2_sequence"]})
            elif scheme == "SNAIL-FISH":
                if p.metadata.get("primer_sequence"):
                    out.append({"Name": f"{p.probe_id}_primer",
                                "Sequence": p.metadata["primer_sequence"]})
                if p.metadata.get("padlock_sequence_5phos"):
                    out.append({"Name": f"{p.probe_id}_padlock",
                                "Sequence": p.metadata["padlock_sequence_5phos"]})
            else:
                out.append({"Name": p.probe_id, "Sequence": p.sequence})
        with open(path, "w", newline="", encoding="utf-8-sig") as fh:
            writer = csv.DictWriter(fh, fieldnames=["Name", "Sequence"])
            writer.writeheader()
            writer.writerows(out)
        self.status_var.set(f"订购表已导出 {len(out)} 行 → {path}")

    def save_report(self):
        if self.result is None:
            messagebox.showinfo("提示", "请先运行一次设计。")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".txt", initialfile="design_report.txt",
            filetypes=(("文本", "*.txt"),),
        )
        if not path:
            return
        Path(path).write_text(self.report_txt.get("1.0", "end"), encoding="utf-8")
        self.status_var.set(f"报告已保存 → {path}")

    # ------------------------------------------------------------------
    # 基因组管理
    # ------------------------------------------------------------------
    def refresh_backgrounds(self):
        for child in self.bg_list_frame.winfo_children():
            child.destroy()
        self.bg_vars.clear()
        registry = load_registry()
        if not registry:
            self.bg_hint.set("尚无已注册基因组；请到「基因组与索引」页注册。")
            return
        for gid in sorted(registry):
            var = tk.BooleanVar(value=False)
            self.bg_vars[gid] = var
            ttk.Checkbutton(
                self.bg_list_frame,
                text=f"{registry[gid].get('organism', gid)}（{gid}）",
                variable=var,
            ).pack(anchor="w")
        self.bg_hint.set("仅\"物种区分\"场景需要勾选；低丰度液培检测保持不勾选。")
        self._refresh_genome_list()

    def _refresh_genome_list(self):
        if not hasattr(self, "genome_list"):
            return  # 右侧页签尚未构建
        self.genome_list.delete(0, "end")
        for gid, value in sorted(load_registry().items()):
            kind = "背景/宿主" if value.get("is_host") else "靶标"
            self.genome_list.insert("end", f"{gid}  {value.get('organism', '')}  [{kind}]")

    def register_genome(self):
        path = filedialog.askopenfilename(
            title="选择基因组 FASTA",
            filetypes=(("FASTA", "*.fa *.fasta *.fna *.txt"), ("所有文件", "*.*")),
        )
        if not path:
            return
        gid = simple_ask_string(self, "索引 ID", "例如 ecoli_K12：")
        if not gid:
            return
        organism = simple_ask_string(self, "显示名称", "例如 E. coli K-12：") or gid
        registry = load_registry()
        if gid in registry and not messagebox.askyesno(
                "覆盖确认", f"索引 {gid} 已存在，是否覆盖？"):
            return
        self.status_var.set(f"正在构建 {gid} 的 bowtie2 索引…")
        self.run_btn.configure(state="disabled")

        def worker():
            try:
                GENOME_DIR.mkdir(parents=True, exist_ok=True)
                INDICES_DIR.mkdir(parents=True, exist_ok=True)
                dst = GENOME_DIR / f"{gid}.fa"
                dst.write_bytes(Path(path).read_bytes())
                prefix = str(INDICES_DIR / gid)
                build_bowtie2_index(str(dst), prefix, threads=2)
                registry = load_registry()
                registry[gid] = {
                    "organism": organism,
                    "index_prefix": prefix,
                    "fasta_path": str(dst),
                    "is_host": True,
                }
                save_registry(registry)
                self.ui_queue.put(("genome_ok", gid))
            except AlignmentError as exc:
                self.ui_queue.put(("error", str(exc)))

        threading.Thread(target=worker, daemon=True).start()

    def delete_genome(self):
        selection = self.genome_list.curselection()
        if not selection:
            messagebox.showinfo("提示", "请先在列表中选中一个基因组。")
            return
        gid = load_registry_sorted()[selection[0]][0]
        if not messagebox.askyesno("删除确认", f"确定删除索引 {gid}？（FASTA 与索引文件保留在磁盘）"):
            return
        registry = load_registry()
        registry.pop(gid, None)
        save_registry(registry)
        self.refresh_backgrounds()


def statistics_mean(values):
    return sum(values) / len(values)


def statistics_sd(values):
    if len(values) < 2:
        return 0.0
    mean = sum(values) / len(values)
    return (sum((x - mean) ** 2 for x in values) / (len(values) - 1)) ** 0.5


def load_registry() -> dict:
    if REGISTRY_PATH.is_file():
        try:
            return json.loads(REGISTRY_PATH.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            return {}
    return {}


def load_registry_sorted():
    return sorted(load_registry().items())


def save_registry(registry: dict) -> None:
    REGISTRY_PATH.write_text(
        json.dumps(registry, ensure_ascii=False, indent=2), encoding="utf-8"
    )


def simple_ask_string(parent, title: str, prompt: str) -> str | None:
    """极简输入对话框（避免引入 tkinter.simpledialog 的额外样式差异）。"""
    dialog = tk.Toplevel(parent)
    dialog.title(title)
    dialog.transient(parent)
    dialog.grab_set()
    result: list[str | None] = [None]

    ttk.Label(dialog, text=prompt, padding=8).pack(anchor="w")
    entry = ttk.Entry(dialog, width=32)
    entry.pack(padx=8, pady=(0, 8))
    entry.focus_set()

    def confirm(_event=None):
        result[0] = entry.get().strip()
        dialog.destroy()

    def cancel(_event=None):
        dialog.destroy()

    btns = ttk.Frame(dialog)
    btns.pack(pady=(0, 8))
    ttk.Button(btns, text="确定", command=confirm).pack(side="left", padx=4)
    ttk.Button(btns, text="取消", command=cancel).pack(side="left", padx=4)
    entry.bind("<Return>", confirm)
    entry.bind("<Escape>", cancel)
    parent.wait_window(dialog)
    return result[0]


def main():
    root = tk.Tk()
    try:
        style = ttk.Style(root)
        style.theme_use("aqua")
    except tk.TclError:
        pass
    ProbeGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
