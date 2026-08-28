# ProbeStudio · FISH 探针设计本地 GUI

ProbeStudio 是一个完全在本机运行的 FISH 探针设计与过滤平台，基于 J 的
探针设计引擎（原 FastAPI/React 方案的核心算法层）完全重构而来：

- **四种方案 · 独立模块**：smFISH、smiFISH、HCR 3.0、SNAIL FISH，各自是
  `probedesign/schemes/` 下的独立设计模块，可扩展新方案。
- **宿主基因组过滤**：候选探针比对到注册的宿主/背景基因组（bowtie2），
  超阈值即淘汰——降低感染/定植模型中的背景信号。
- **物理模型**：SantaLucia 1998 最近邻 Tm（单价盐 + Mg²⁺/dNTP + 甲酰胺
  校正）、primer3 发卡、Sugimoto 1995 RNA/DNA Gibbs 自由能（HCR）。
- **界面**：侧边栏参数全部带解释；结果含指标卡、漏斗图、覆盖图、分布图、
  单条详情与订购表导出；不做"合格/不合格"自动判定。

## 环境

- conda 环境：`/opt/anaconda3/envs/Probe`（Python 3.11）
- 依赖：streamlit、pandas、numpy、biopython、primer3-py、plotly、bowtie2

安装（在 Probe 环境内）：

```bash
pip install -e .
conda install -n Probe -c bioconda -c conda-forge bowtie2
```

## 启动

双击 `launch.command`，或在终端：

```bash
/opt/anaconda3/envs/Probe/bin/python -m streamlit run app.py
```

浏览器会自动打开本地界面（仅监听 127.0.0.1，序列不出机）。

## 使用流程

1. **注册基因组索引**（一次性）：展开"基因组与索引管理"，上传目标/宿主
   FASTA 构建 bowtie2 索引；或把现成的索引目录写入 `genome_registry.json`。
2. **侧边栏**：粘贴靶序列 FASTA（或选择已注册靶标索引）→ 选择设计方案
   → 调整方案参数（每项带解释）→ 选择宿主/背景基因组 → 开始设计。
3. **结果区**：指标卡（候选数/通过数/覆盖率/Tm 均一性）→ 漏斗图 →
   结果表（含方案专属列）→ 覆盖图/分布图 → 单条详情 → 导出
   （CSV / IDT 兼容订购表 / 参数汇总 JSON）。

## 方案说明

| 方案 | 产物 | 关键过滤 |
|------|------|----------|
| smFISH | 18–24 nt 反义寡核苷酸 | Tm/GC 窗口、发卡 Tm、同聚碱基、特异性 |
| smiFISH | smFISH + 共享 readout 延伸段 | 同上 + 完整序列组装 |
| HCR 3.0 | 两条 25-mer 半探针（带分裂 initiator） | tile GC、Gibbs 窗口、半探针 dTm |
| SNAIL FISH | primer + 5′-磷酸化 padlock（含 UGI 条码） | 双臂独立 GC/重复/发卡 dG、三级特异性 |

## 代码结构

```
src/probedesign/
├── models.py      # 数据模型（参数/探针/结果）
├── config.py      # SantaLucia NN 表与默认条件
├── utils.py       # Tm、GC、反向互补等序列工具
├── mining.py      # 候选枚举
├── filters.py     # 热力学过滤（廉价过滤优先，发卡最后）
├── alignment.py   # bowtie2 封装与 SAM 计数
├── scoring.py     # 特异性过滤与排序
├── selection.py   # 间距与降采样
├── report.py      # 输出表
├── pipeline.py    # run_design 入口
└── schemes/       # smfish / smifish / hcr3 / snail 独立模块
app.py             # Streamlit GUI
```

## 与旧版（FastAPI/React）相比的引擎修正

- Tm 最近邻表修正为完整 SantaLucia 1998（原表混入错误二核苷酸值、缺
  AA/TT/GG/CC、末端 AT 罚分符号反了）；与 primer3 独立实现偏差 <3 °C。
- bowtie2 比对计数纳入 secondary 记录（旧解析只数 primary，每条探针
  最多计 1 次，重复/宿主过滤实际失效）。
- 探针间距改为真正的区间重叠判断（旧版只比较起点距离，会放过重叠探针）。
- HCR 3.0：半探针拆分与 dTm 过滤移到选点之前（旧版可能选中事后被淘汰
  的探针）；tile Tm 窗口默认关闭（52-mer 必然超出 smFISH Tm 窗口，
  协议以 dTm + Gibbs 为准）。
- SNAIL：primer/padlock 的结合臂改为反义链（旧版用了正链序列，无法
  与靶 RNA 杂交）；补算候选 Tm/GC。
- 热力学过滤按"廉价优先"排序，primer3 发卡只在存活候选上计算（大幅提速）。
- bowtie2 缺失/失败时报错明确（旧版吞掉 stderr）。

## 测试

```bash
PYTHONPATH=src /opt/anaconda3/envs/Probe/bin/python -m pytest tests/ -q
```

`data/smoke/` 内含 500 nt 合成靶标 + 宿主基因组及其 bowtie2 索引，可直接
在界面中注册 `data/smoke/host_idx`、`data/smoke/target_idx` 做端到端体验。

## 许可证

MIT
