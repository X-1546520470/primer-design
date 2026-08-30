# MycoPrimerV2 · 分枝杆菌 FISH 探针设计 GUI

MycoPrimerV2（原 ProbeStudio / primer-design v1）是一个完全在本机运行的
分枝杆菌 FISH 探针设计与过滤平台，V2 按三轮批量测试的结论对四种 FISH
方案的脚本做了独立优化，并整合进同一个 GUI：

- **四种方案 · 独立模块**：smFISH、smiFISH、HCR 3.0、SNAIL FISH，各自的
  设计脚本独立实现（`src/mycoprimer/schemes/`），注册表分发、可扩展。
- **设计目标预设**（V2 新增）：低丰度转录本检测（液培）与物种区分探针
  两套经验配置一键套用——来源见仓库内三份批量测试报告。
- **分枝杆菌调优的 HCR 窗口**（V2 新增）：针对 65–67% GC 基因组将默认
  窗口调整为 GC 40–65、dTm≤8、Gibbs −75~−45（官方中低 GC 窗口在
  分枝杆菌上 48 基因测试仅 1.4 对/基因，调参后 3.9 对/基因）。
- **SNAIL 订购便捷**（V2 新增）：padlock 自动生成带 /5Phos/ 的订购变体。
- **smiFISH 校验**（V2 新增）：readout/linker 含非 ACGT 字符直接报错。
- **宿主基因组过滤**：候选探针比对到注册的背景基因组（bowtie2），按需开关。
- **本地运行**：仅监听 127.0.0.1，序列不出机；所有指标为描述性参考，
  不做"合格/不合格"自动判定。

V1 → V2 的引擎修正（含 5 个批量测试发现的 bug）见下文"引擎修正"。

## 环境

- conda 环境：`/opt/anaconda3/envs/Probe`（Python 3.11）
- 依赖：streamlit、pandas、numpy、biopython、primer3-py、plotly、bowtie2

安装（在 Probe 环境内）：

```bash
pip install -e .
conda install -n Probe -c bioconda -c conda-forge bowtie2
```

## 启动

### 在新的电脑上（首次）

1. `git clone -b mycoprimer-v2 https://github.com/X-1546520470/primer-design.git`
   （或下载 zip 解压）
2. 双击 **`setup.command`** —— 自动完成：Python ≥3.10 环境创建、
   primer3-py/biopython 等依赖安装、bowtie2 安装（conda 或 Homebrew）。
3. 双击 **`run_gui.command`** 打开桌面版探针设计界面。

要求：macOS；Python ≥3.10（setup 会自动寻找）；bowtie2（setup 自动装，
装不上时按提示手动装）。

### 日常启动

- **桌面版（推荐，Tkinter 原生窗口）**：双击 `run_gui.command`
- **网页版（Streamlit，功能相同）**：双击 `launch.command`
- 命令行批量：`mycoprimer-batch 基因.fasta --scheme smFISH --out-dir 结果目录`

界面布局：① 靶序列输入 → 左侧方案与参数（设计目标预设一键套用）→
右侧页签（结果表 / 单条详情 / 设计报告 / 导出 / 基因组与索引 / 使用说明）。
多记录 FASTA（多条 `>基因名`）会自动逐基因批量设计。

## 启动（旧）

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
src/mycoprimer/
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

📖 详细操作指引见 [使用说明.md](使用说明.md)（含 MTB/BCG/MSM 基因组数据说明与设计流程）。
