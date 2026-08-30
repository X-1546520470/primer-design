# 丙酰辅酶A（propionyl-CoA）代谢基因 · 四套探针设计工作流批量测试报告

- **测试日期**：2026-08-29
- **测试对象**：从 MycoBrowser / NCBI RefSeq 注释中选取的 **40 个丙酰辅酶A
  代谢相关基因**（M. tuberculosis H37Rv 16 + M. smegmatis mc²155 15 +
  M. bovis BCG Pasteur 1173P2 9），要求 ≥30，实际 40。
- **被测工作流**：**smFISH、smiFISH、HCR 3.0、SNAIL-FISH** 四套（`src/mycoprimer/schemes/`）。
- **测试规模**：4 方案 × 40 基因 = **160 次完整设计运行**，全部端到端
  （候选枚举 → 热力学过滤 → bowtie2 靶标/背景比对 → 打分 → 间距选点 → 产物组装/自检）。
- **测试脚本**：[batch_propionyl_matrix_test.py](batch_propionyl_matrix_test.py)
  （**只测试、不修改引擎**）；明细数据见
  [test_data/propionyl_matrix_results/propionyl_matrix_summary.csv](test_data/propionyl_matrix_results/propionyl_matrix_summary.csv)。
- **重要说明**：本轮发现的所有脚本/引擎问题**均未修正**，按要求统一记录在
  第 5 节"发现的系统性问题（待修复）"，供后续处理。

---

## 1. 总体结论

| 方案 | 运行 | 崩溃/异常 | 结构自检 | 总探针 | 平均/基因 | 零产出基因 | 结论 |
|------|------|-----------|----------|--------|-----------|------------|------|
| **smFISH** | 40 | 0 | 40/40 通过 | 291 | 7.3 | 13 | ✅ 流程稳定；零产出集中在 MTB/BCG（跨物种背景过滤，预期行为） |
| **smiFISH** | 40 | 0 | 40/40 通过 | 291 | 7.3 | 13 | ✅ 与 smFISH 完全一致（复用同一引擎），readout 拼接正确 |
| **HCR 3.0** | 40 | 0 | 40/40 通过 | 60 | 1.5 | 25 | ⚠️ 流程可用；高 GC 下 tile GC/Gibbs 窗口是硬瓶颈，产出远低于目标 |
| **SNAIL-FISH** | 40 | 0 | 40/40 通过 | 104 | 2.6 | 25 | ⚠️ 流程可用；双臂 GC/发卡过滤产率低；**第三级 padlock 特异性检查在 N 占位时失效（问题 P3）** |

**一句话总结**：四套工作流在 160 次运行中**零崩溃、零异常、产物结构自检
全部通过**，管线（含 bowtie2 比对、打分、选点、序列组装）端到端可用。
零产出主要来自两个非崩溃性原因：① 用另外两个物种做背景过滤时，**MTB 与
BCG 同源度 >99.9%，保守代谢基因的探针几乎被背景过滤全部淘汰**（工具的
预期行为）；② HCR/SNAIL 在 65–67% GC 的分枝杆菌基因组上热力学窗口偏严，
MSM 上也仅 4.0 / 6.9 条每基因。另发现 4 个需要修复的系统性问题（第 5 节）。

---

## 2. 基因选取与数据来源

### 2.1 来源

| 物种 / 株 | 注释来源 | 基因组 FASTA / bowtie2 索引 | 命中丙酰辅酶A基因 | 实际取用 |
|-----------|----------|------------------------------|--------------------|----------|
| M. tuberculosis H37Rv | **MycoBrowser** GFF（`mycobrowser_h37rv.gff`） | `genomes/mtb_h37rv.fna` / `indices/mtb_h37rv` | 17 | 16 |
| M. smegmatis mc²155 | **MycoBrowser** GFF（`mycobrowser_smegmatis.gff`） | `genomes/msm_mc2155.fna` / `indices/msm_mc2155` | 15 | 15 |
| M. bovis BCG Pasteur 1173P2 | **NCBI RefSeq** GFF（`GCF_000009445.1_...gff`；MycoBrowser 不收录 BCG，按 H37Rv 同源注释取） | `genomes/bcg_pasteur.fna` / `indices/bcg_pasteur` | 9 | 9 |
| **合计** | | | **41** | **40** |

### 2.2 筛选标准（满足其一，见 [batch_smfish_test.py](batch_smfish_test.py)）

- a) 功能注释含 `propionyl / methylcitrate / methylisocitrate / methylmalonyl`
  （丙酰辅酶A直接代谢酶）；
- b) 基因名属 `prpB / prpC / prpD / prpR`（甲基柠檬酸循环）；
- c) 基因名属 `icl1 / icl2 / aceAa / aceAb / aceA`（乙醛酸支路，MTB 中
  丙酰辅酶A脱毒的替代通路）。

覆盖的功能包括：甲基柠檬酸循环（prpC/prpD/prpB）、甲基丙二酰辅酶A变位酶
（mutA/mutB/scpA）、甲基丙二酰辅酶A差向酶（mce）、乙酰/丙酰辅酶A羧化酶
各亚基（accA1/A2/A3、accD1/D2/D5/D6、accE5）、异柠檬酸裂解酶（icl1）、
α-甲基酰基辅酶A消旋酶（mcr）、MeaB 等。基因 FASTA 见
[test_data/batch_results/gene_fastas/](test_data/batch_results/gene_fastas/)，
清单与产物见 [test_data/batch_results/batch_summary.csv](test_data/batch_results/batch_summary.csv)。

### 2.3 参数与背景设置（本轮）

- 各方案参数：smFISH/smiFISH/HCR3/SNAIL 均 `desired_probe_count=20`、
  `min_gap=2`；smiFISH 用占位 readout `ACGTCGACTATCGAT`（3′ 端）；
  HCR3 用通道 B1、分枝杆菌调优窗口（GC 40–65、dTm≤8、Gibbs −75~−45）；
  SNAIL 用默认 20 nt 双臂、UGI 以 22 个 N 占位。
- **背景基因组（交叉反应过滤）**：设计某物种探针时用另外两个物种做背景——
  MTB→{BCG,MSM}，BCG→{MTB,MSM}，MSM→{MTB,BCG}。
- bowtie2：`--very-sensitive-local`，`--score-min C,36,0`，`-k 100`，4 线程；
  计数含 secondary 比对。

---

## 3. 总体结果（按方案 × 物种）

### 3.1 漏斗汇总（候选 → 热力学存活 → 特异性存活 → 最终，各格为该组均值）

| 方案 | 物种 | 候选 | 热力学存活 | 特异性存活 | 背景淘汰 | 靶重复淘汰 | 最终探针 |
|------|------|------|------------|------------|----------|------------|----------|
| smFISH | MTB | 9 998 | 1 165 | 11 | 1 154 | 0 | **0.6** |
| smFISH | BCG | 10 768 | 1 147 | 5 | 1 141 | 0 | **0.4** |
| smFISH | MSM | 10 566 | 826 | 786 | 40 | 0 | **18.5** |
| smiFISH | MTB/BCG/MSM | 同 smFISH | 同 smFISH | 同 smFISH | 同 smFISH | 0 | **0.6 / 0.4 / 18.5** |
| HCR3 | MTB | 1 397 | 72 | 0 | 72 | 0 | **0.0** |
| HCR3 | BCG | 1 507 | 58 | 0 | 58 | 0 | **0.0** |
| HCR3 | MSM | 1 478 | 29 | 24 | 5 | 0 | **4.0** |
| SNAIL | MTB | 1 408 | 130 | 0 | 130 | 0 | **0.0** |
| SNAIL | BCG | 1 518 | 110 | 0 | 99 | 0 | **0.0** |
| SNAIL | MSM | 1 489 | 88 | 74 | 13 | 0 | **6.9** |

### 3.2 背景淘汰比例（热力学存活候选中被跨物种背景淘汰的比例）

| 方案 | MTB | BCG | MSM |
|------|-----|-----|-----|
| smFISH/smiFISH | **99.1%** | **99.5%** | 4.9% |
| HCR3 | **100%** | **100%** | 18.0% |
| SNAIL | **100%** | **100%** | 15.2% |

### 3.3 覆盖率与 Tm（仅有产出的基因）

| 方案 | 有产出基因 | 覆盖率均值 | 覆盖率中位 | 报告 Tm 均值范围 (°C) |
|------|------------|------------|------------|------------------------|
| smFISH | 27 | 14.1% | 17.0% | 63.1–69.7（18–24 nt 探针 Tm，可横向比较） |
| smiFISH | 27 | 14.1% | 17.0% | 63.1–69.7（同 smFISH） |
| HCR3 | 15 | 14.7% | 15.5% | 87.3–90.2（**52 nt tile 整体 Tm**，非半探针 Tm，不可与 smFISH 横比，见问题 P6） |
| SNAIL | 15 | 19.2% | 18.6% | 83.7–86.0（**41 nt 双臂 cassette Tm**，非 20 nt 臂 Tm，不可横比，见问题 P6） |

---

## 4. 关键现象与归因

### 4.1 MTB/BCG 零产出 = 跨物种背景过滤的预期行为（非缺陷）

丙酰辅酶A代谢是**高度保守的核心代谢**。MTB 与 BCG 基因组同源度 >99.9%，
针对 MTB/BCG 保守基因设计的探针几乎必然在另一物种基因组上命中。背景过滤
（`max_host_hits=0` 零容忍）因此把 99–100% 的热力学存活候选淘汰：

- smFISH/smiFISH：MTB 零产出 8/16、BCG 5/9；MSM 0/15。
- HCR3/SNAIL：MTB 16/16、BCG 9/9 全部零产出；MSM 0/15。

这说明这些保守基因**不适合用作"MTB vs BCG 物种区分"探针**——这正是第一
轮物种区分测试的结论，本轮在四套方案上一致复现。若实验目的是**纯培养
独立检测**（不做物种区分），应像低丰度测试那样关闭跨物种背景（hosts=[]），
此时 MTB/BCG 同样可产出（参见低丰度报告）。**这是配置/生物学问题，不是
脚本 bug。**

### 4.2 HCR 3.0 高 GC 瓶颈：tile GC 窗口主导（MSM 15 基因漏斗拆解）

对 MSM 15 基因的 22 176 个 tile 候选按**首要失败原因**归类：

| 淘汰阶段 | 数量 | 占比 |
|----------|------|------|
| tile GC 超窗（40–65%） | 17 290 | **78.0%** |
| 半探针发卡 Tm 超标 | 3 267 | 14.7% |
| Gibbs 自由能超窗 | 1 029 | 4.6% |
| 未入选（间距/降采样） | 370 | 1.7% |
| 同聚碱基 | 130 | 0.6% |
| dTm 超标 | 27 | 0.1% |
| **最终保留** | **63** | **0.3%** |

即便已用分枝杆菌调优窗口（GC 40–65，较官方 45–55 放宽），仍有 78% 候选
因 tile GC 出局。这与既往 48 基因测试结论一致（[三方案可行性对比测试报告.md](三方案可行性对比测试报告.md)）：
**高 GC 序列组成是硬约束，不是工具问题**。HCR 每对探针启动一条放大器链、
信号指数放大，MSM 上 4.0 对/基因通常仍可检出；要提高对数需放宽窗口或扩展
靶区（含 UTR）。

### 4.3 SNAIL 高 GC 瓶颈：双臂 GC 窗口主导（MSM 15 基因漏斗拆解）

对 MSM 15 基因的 22 341 个双臂候选按首要失败原因归类：

| 淘汰阶段 | 数量 | 占比 |
|----------|------|------|
| 单臂 GC 超窗（40–63%） | 20 677 | **92.6%** |
| 未入选（间距/降采样） | 1 202 | 5.4% |
| 臂重复基序（AAAA/CCCC/GGGG/TTTT） | 350 | 1.6% |
| **最终保留** | **112** | **0.5%** |

SNAIL 要求**两条 20 nt 臂同时**落在 40–63% GC 窗口，在 ~65% GC 基因组上
概率极低，92.6% 候选被臂 GC 淘汰。这是方案本身的序列约束所致。

### 4.4 smFISH 与 smiFISH 结果完全一致

smiFISH 复用 `design_smfish` 的全部流程，仅在末尾拼接 readout。本轮两者
候选数、漏斗各阶段、最终 291 条探针、覆盖率、Tm **逐基因完全相同**，
smiFISH 的 `full_sequence = 探针 + TTT + readout` 结构自检 40/40 通过，
readout/linker 非 ACGT 字符校验逻辑在位（未触发，因占位序列合法）。

### 4.5 靶标自身重复序列过滤本轮未触发

所有 160 次运行 `target_rejected`（靶基因组命中 >10）均为 **0**。丙酰辅酶A
基因多为单拷贝核心酶，符合预期；该 QC 通道在本基因集上未被压力测试到。

---

## 5. 发现的系统性问题（**按要求暂不修正，仅记录**）

> 以下问题均在本轮测试中观察/验证。除 P1 影响所有运行外，其余不影响
> 本轮"管线能跑通"的结论，但建议在后续版本处理。

### P1 · FASTA 文件句柄未关闭（ResourceWarning，全部 160 次运行均触发）

- **现象**：每次设计运行都产生 1 条 `ResourceWarning: unclosed file <_io.TextIOWrapper ...fasta mode='rt'>`。
- **定位**：
  - HCR3/SNAIL 走 [schemes/common.py `load_first_target`](src/mycoprimer/schemes/common.py#L83-L89)：
    `records = list(SeqIO.parse(target_fasta, "fasta"))` 直接传文件**路径**，
    Biopython 内部打开的句柄未显式关闭（警告定位在 `common.py:86`）。
  - smFISH/smiFISH 走 [mining.py `load_fasta`](src/mycoprimer/mining.py#L22-L24)：
    同样 `list(SeqIO.parse(path, "fasta"))` 传路径。
- **影响**：低。句柄依赖 GC 回收关闭，单次设计无碍；但批量/长会话（GUI
  多次设计）下会累积打开句柄，属于资源泄漏隐患。
- **建议（不在本轮改）**：改为 `with open(path) as fh: list(SeqIO.parse(fh, "fasta"))`，
  或用 `with` 上下文管理。

### P2 · SNAIL 选点间距存在"区间长度双重计数"（latent，本轮数据上未压低产量）

- **现象/定位**：[schemes/snail.py `design_snail`](src/mycoprimer/schemes/snail.py#L246-L253)
  计算 `min_span = 2*arm_len + spacer + min_gap`（默认 41+2=43），再把它
  作为 `min_gap` 传给 [selection.py `select_non_overlapping`](src/mycoprimer/selection.py#L21-L23)。
  而 `_intervals_too_close` 判定的是 `a.start < b.stop + min_gap`，
  **探针自身区间长度（cassette 41 nt）已经包含在 `b.stop−b.start` 里**。
  于是实际要求的相邻起点间距 ≈ 区间长 41 + min_gap 43 = **84 nt**，是预期
  "相邻结合区留 2 nt 间隔（起点间距 43 nt）"的约两倍。
- **本轮影响**：在丙酰辅酶A基因上**实测无差异**——存活候选本就稀疏聚集，
  用 `min_gap=2` 与 `min_gap=43` 重选，MSM 15 基因逐条产量完全相同
  （11/8/10/6/11/7/6/1/11/3/8/8/7/7/8）。因此属**潜在逻辑问题**，在候选
  密集的长转录本/低 GC 靶标上可能使 SNAIL 产量减半。
- **建议（不在本轮改）**：SNAIL 选点应直接传 `min_gap=params.min_gap`
  （区间长度已由选点函数内含），或改用"起点间距"语义的独立判定。

### P3 · SNAIL 第三级 padlock 特异性检查在 UGI 用 N 占位时**完全失效**（重要）

- **现象**：40 个基因的 SNAIL 运行中，**从未出现任何 `primer_*` /
  `padlock_*` 级别的特异性淘汰**。
- **根因**：[schemes/snail.py `_assemble_oligos`](src/mycoprimer/schemes/snail.py#L135)
  在未提供真实 UGI 条码时用 22 个连续 N 占位
  （`ugi = params.snail_ugi_sequence or "N"*22`）。随后
  [`_check_component_specificity`](src/mycoprimer/schemes/snail.py#L153-L202)
  把含 22 N 的 padlock 提交 bowtie2 比对。bowtie2 将 N 视为不可匹配
  （ambiguous，不计分且禁止落在 seed 区），22 个连续 N 横跨序列中部，
  导致 padlock **永远比对不上、命中数恒为 0**，第三级 padlock 检查因此
  从不淘汰任何候选。
- **直接验证**：同一条 SNAIL 产物，padlock（含 22 N）比对 MSM 索引命中
  **0**；primer（无 N）命中 **1**。
- **影响**：这是**测试/占位场景下的假阴性**——第三级特异性给出"全部
  通过"的假象，实际并未检查 padlock。使用**真实条码（无 N）时该检查恢复
  有效**，故不影响最终订购序列的安全性，但会让工作流测试/QC 结果过于乐观。
- **建议（不在本轮改）**：比对前把 padlock 中的 N 占位段剔除/截断，或对
  含 N 组件跳过 bowtie2 并显式标注"占位未检"，避免静默失效。

### P4 · 失败原因分类中 HCR 的 Gibbs 淘汰未被归入"热力学"（统计口径，非引擎崩溃）

- **现象**：本轮测试脚本的漏斗阶段分类（以及 [batch_smfish_test.py](batch_smfish_test.py)
  的 `failure_category`）用关键词 `GC=/Tm=/homopolymer/hairpin` 识别热力学
  淘汰。HCR tile 的 **Gibbs 自由能**失败原因串是 `Gibbs=... outside [...]`，
  不含上述关键词，会被归入 "other" 而非 "thermo"。
- **影响**：仅影响**测试/报告脚本的漏斗计数口径**，不影响引擎本身的过滤
  行为（引擎正确淘汰了这些候选）。本轮报告的 HCR 漏斗数字已用独立脚本按
  完整关键词（含 Gibbs/dTm/arm）重新归类（见 4.2/4.3）。
- **建议（不在本轮改）**：在共用的 `failure_category` 关键词中补 `Gibbs`、
  `dTm`、`arm1/arm2`，让四方案的漏斗统计口径统一。

### P5 · smFISH/smiFISH 对已热力学淘汰的候选仍提交比对（效率不一致，非错误）

- **现象**：[schemes/smfish.py](src/mycoprimer/schemes/smfish.py#L50-L64)
  把**全部候选**（含热力学已淘汰者）都提交 bowtie2；而 HCR3/SNAIL 走
  [common.py `apply_target_alignment`/`apply_host_alignment`](src/mycoprimer/schemes/common.py#L101-L158)，
  只对 `p.passed` 的存活候选比对。
- **影响**：结果正确（smFISH 后续 `apply_specificity_filters` 只处理存活
  候选），但 smFISH/smiFISH 比对了大量本会被热力学淘汰的序列，bowtie2
  调用量偏大（本轮 smFISH 每基因约 1 万候选全部比对）。属实现风格/效率
  不一致，不影响正确性。
- **建议（不在本轮改）**：统一为"只比对存活候选"，减少比对规模。

### P6 · 报告 Tm 字段在不同方案间语义不同，易被误读（展示层）

- **现象**：结果表/汇总里的 `tm_mean`，对 smFISH/smiFISH 是 18–24 nt 探针
  Tm（63–70 °C）；对 HCR3 是 **52 nt tile 整体 Tm**（87–90 °C）；对 SNAIL
  是 **41 nt 双臂 cassette Tm**（83–86 °C）。三者尺度不同、**不可横向比较**。
- **影响**：GUI/报告若把四方案 Tm 画在同一轴或直接对比，会得出"HCR/SNAIL
  Tm 异常高"的错误印象。HCR 真正有意义的是两条 25-mer 半探针 Tm 与 dTm；
  SNAIL 是两条 20 nt 臂的 Tm。
- **建议（不在本轮改）**：结果模型按方案暴露"实际合成序列的 Tm"
  （HCR 半探针 Tm、SNAIL 臂 Tm），并在 GUI/报告中标注 Tm 对应的序列范围。

---

## 6. 复现方法

```bash
# 在项目目录、Probe conda 环境下
/opt/anaconda3/envs/Probe/bin/python batch_propionyl_matrix_test.py
# 汇总输出：test_data/propionyl_matrix_results/propionyl_matrix_summary.csv
```

- 基因选取/序列提取逻辑复用第一轮 [batch_smfish_test.py](batch_smfish_test.py)
  （MycoBrowser GFF + NCBI RefSeq GFF），可独立重跑复现 40 基因。
- 本轮调查脚本（漏斗拆解、间距对照、padlock-N 验证）为只读分析，未改动
  引擎；四套方案模块与管线代码在测试前后保持不变。

---

## 7. 建议的后续动作（按优先级）

1. **P3（SNAIL padlock N 占位致第三级特异性静默失效）**——影响 QC 可信度，
   建议优先修：比对前处理 N 占位或显式标注"占位未检"。
2. **P1（FASTA 句柄未关闭）**——批量/GUI 长会话的资源泄漏，修复成本低。
3. **P2（SNAIL 间距双重计数）**——当前数据无影响，但建议修正以避免低 GC/
   长靶标上产量被腰斩。
4. **P4/P6（统计口径与 Tm 展示）**——统一漏斗分类关键词、按方案标注 Tm
   语义，提升报告/GUI 可读性。
5. **P5（smFISH 比对效率）**——统一只比对存活候选。
6. 若实验目标是**纯培养检测而非物种区分**，丙酰辅酶A这类保守基因应
   **关闭跨物种背景过滤**（hosts=[]），否则 MTB/BCG 上会被背景过滤清零；
   HCR/SNAIL 在高 GC 靶标上建议进一步放宽窗口或扩展靶区以提高对数。
