# 丙酰辅酶A（propionyl-CoA）代谢基因 · 低丰度转录检测 · 四方案批量测试报告

- **测试日期**：2026-08-29
- **实验目的**：在**液体培养、低菌量**条件下，**独立检测**几种分枝杆菌
  （M. tuberculosis H37Rv、M. bovis BCG Pasteur、M. smegmatis mc²155）
  **低丰度靶基因的转录表达**。各菌种**独立培养、独立制片检测**，样本中
  只含该菌种自身 RNA。
- **关键配置变更**（相对上一轮"物种区分"测试）：
  - **关闭跨物种背景过滤**（`host_genomes=[]`）——不做物种区分，也不应
    用另一物种基因组淘汰探针（否则保守代谢基因会被误杀）；
  - **保留靶标自身基因组的重复序列 QC**（`max_target_hits=10`）——这是
    探针质量的内在检查，与跨物种过滤无关；
  - 低丰度检测提高单转录本探针堆叠：smFISH/smiFISH 目标 **48** 条
    （经典 smFISH 高堆叠），HCR3 **20** 对、SNAIL **20** 对（两者信号
    经扩增放大，无需 48 条）。
- **基因集**：MycoBrowser + NCBI RefSeq 注释中**重新选取**的 **40 个
  丙酰辅酶A代谢相关基因**（MTB 16 + MSM 15 + BCG 9，≥30）。
- **测试规模**：4 方案 × 40 基因 = **160 次端到端设计**。
- **脚本/数据**：[batch_propionyl_lowabundance_test.py](batch_propionyl_lowabundance_test.py)；
  明细 [propionyl_lowabundance_summary.csv](test_data/propionyl_lowabundance_results/propionyl_lowabundance_summary.csv)；
  对比表 [丙酰辅酶A低丰度四方案对比明细表.md](丙酰辅酶A低丰度四方案对比明细表.md)。
- 本轮仍**只测试、不修改引擎**。

---

## 1. 总体结论

| 方案 | 运行 | 崩溃 | 零产出 | 结构自检 | 总探针 | 均值/基因 | 范围 | 覆盖率均值 |
|------|------|------|--------|----------|--------|-----------|------|-----------|
| **smFISH** | 40 | 0 | **0** | 40/40 通过 | 1397 | 34.9 | 7–48 | 43.2% |
| **smiFISH** | 40 | 0 | **0** | 40/40 通过 | 1397 | 34.9 | 7–48 | 43.2% |
| **HCR 3.0** | 40 | 0 | **0** | 40/40 通过 | 221 | 5.5 | 1–11 | 20.7% |
| **SNAIL-FISH** | 40 | 0 | **0** | 40/40 通过 | 308 | 7.7 | 1–11 | 22.1% |

**核心结论**：关闭跨物种背景过滤后，**40 个基因在四套方案下全部成功
产出探针，160 次运行零崩溃、零零产出、零结构问题**。这验证了上一轮
MTB/BCG 的零产出确实是"跨物种背景过滤"这一配置造成的，而非基因本身或
引擎缺陷。

- **smFISH / smiFISH**：最适合低丰度场景。38/40 基因达到 ≥20 条（仅 2 个
  短基因 Rv1489A 231nt=7 条、MSMEG_4880 399nt=8 条因长度不足），多数
  31–48 条，覆盖率 31–61%。高堆叠可在低菌量下提供足够荧光点数。
- **HCR 3.0 / SNAIL**：受高 GC 热力学窗口限制，对数普遍偏少（均值
  5.5 / 7.7 对）。但二者信号经杂交链反应 / 滚环扩增**指数级放大**，单转录本
  仅需少量探针即可检出，且非零产出覆盖全部基因；可作为"低探针数 + 高信号"
  的备选方案，短基因/特高 GC 基因上对数偏少需关注。

---

## 2. 基因选取（重新选取，与第一轮同标准）

| 物种 / 株 | 注释来源 | 命中 | 取用 |
|-----------|----------|------|------|
| M. tuberculosis H37Rv | MycoBrowser GFF | 17 | 16 |
| M. smegmatis mc²155 | MycoBrowser GFF | 15 | 15 |
| M. bovis BCG Pasteur 1173P2 | NCBI RefSeq GFF | 9 | 9 |
| **合计** | | **41** | **40** |

筛选标准：注释含 `propionyl / methylcitrate / methylisocitrate /
methylmalonyl`，或基因名属 `prpB/C/D/R`、`icl1/icl2/aceA*`。覆盖甲基柠檬酸
循环（prpC/prpD/prpB）、甲基丙二酰辅酶A变位酶（mutA/mutB/scpA）、差向酶
（mce）、乙酰/丙酰辅酶A羧化酶亚基（accA/D/E 系列）、异柠檬酸裂解酶
（icl1）、消旋酶（mcr）、MeaB 等。序列见
[test_data/propionyl_lowabundance_results/gene_fastas/](test_data/propionyl_lowabundance_results/gene_fastas/)。

---

## 3. 关键结果

### 3.1 全部基因零零产出

关闭跨物种背景后，连最保守的 icl1、prpC、mutA、mcr、accA1 等上一轮
"四方案全零"的基因，本轮**全部有产出**。例如：

| 基因 | 物种 | 长度 | smFISH | HCR3 | SNAIL |
|------|------|-----:|-------:|-----:|------:|
| icl1 (Rv0467) | MTB | 1287 | 37 | 4 | 7 |
| prpC (Rv1131) | MTB | 1182 | 31 | 3 | 9 |
| mutA (Rv1492) | MTB | 1848 | 33 | 4 | 6 |
| accA1 (Rv2501c) | MTB | 1965 | 43 | 11 | 10 |

靶标自身重复序列淘汰（`target_rejected`）160 次运行**全部为 0**，说明这些
核心代谢基因在各自基因组内均为单拷贝/低拷贝，探针特异性无内在问题。

### 3.2 低产基因（低于低丰度推荐阈值，非零产出）

实用阈值（单转录本信号量经验值）：smFISH/smiFISH ≥20 条、HCR3 ≥8 对、
SNAIL ≥8 对。低于阈值的基因：

- **smFISH/smiFISH（2 个，均为短基因）**：
  - MTB **Rv1489A**（231 nt）= 7 条；MSM **MSMEG_4880**（399 nt）= 8 条。
  - 短转录本物理上放不下更多 18–24 nt 探针，属长度限制，非工具问题。
- **HCR3（多数基因 <8 对）**：高 GC 下 tile GC/Gibbs 窗口淘汰了 ~98%
  候选（详见 3.3），均值仅 5.5 对；≥8 对的多为较长且 GC 适中的基因
  （prpD 9、accA1 11、accD5 9、accD1 9 等）。
- **SNAIL（少数 <8 对）**：均值 7.7 对，短基因 Rv1489A / MSMEG_4880
  仅 1 对，其余多为 5–11 对。

> 完整逐基因数字与 ⚠️ 标注见
> [丙酰辅酶A低丰度四方案对比明细表.md](丙酰辅酶A低丰度四方案对比明细表.md)。

### 3.3 高 GC 仍是 HCR/SNAIL 的产量瓶颈（但非阻断）

在无背景过滤、候选全部可用的情况下，HCR/SNAIL 产量仍显著低于 smFISH，
瓶颈是分枝杆菌 65–67% GC 基因组下的**热力学窗口**：

- HCR3：tile GC 窗口（40–65%）+ Gibbs 自由能窗口 + 半探针发卡/dTm 过滤，
  仅 ~0.3–0.5% 候选存活（与既往报告一致）。
- SNAIL：要求两条 20 nt 臂**同时**落在 40–63% GC 窗口，~92% 候选被臂 GC
  淘汰。

这是方案的序列组成约束，不是脚本错误。由于 HCR/SNAIL 信号经扩增放大，
实际检出能力不完全取决于探针对数；但若需更多对数，可：放宽 GC/Gibbs
窗口（需实验验证非特异结合）、扩展靶区（含 UTR）、或对极短基因改用
smFISH/smiFISH。

### 3.4 smFISH 与 smiFISH 结果完全一致

两者候选、漏斗、最终 1397 条、覆盖率、Tm 逐基因完全相同；smiFISH 的
`探针+TTT+readout` 拼接结构自检 40/40 通过。smiFISH 用共享 readout，
合成成本更低，是低丰度多样本检测的高性价比选择。

---

## 4. 针对实验目的的方案建议

1. **首选 smFISH / smiFISH**：低菌量、低丰度转录检测依赖单转录本上的荧光
   点数，高堆叠（多数基因 31–48 条）信噪比最好；smiFISH 共享 readout 可
   降低多基因/多通道的合成与二级探针成本。两者在本基因集上表现完全等价。
2. **HCR 3.0 / SNAIL 作为备选/补充**：信号扩增带来高亮度，适合转录本极短
   或荧光显微灵敏度受限的情形；但高 GC 下对数偏少，建议对目标基因先看
   明细表中的对数，<8 对的基因谨慎单用，或与 smFISH 组合。
3. **两个极短基因**（Rv1489A 231nt、MSMEG_4880 399nt）：任何方案探针数
   都受限，建议纳入相邻 UTR/共转录区延长靶区，或接受较少探针 + 扩增方案。
4. **订购前必做**（本测试用占位序列）：
   - smiFISH 把占位 readout `ACGTCGACTATCGAT` 替换为实际 LNA 二级探针
     互补序列；
   - HCR3 选定通道（本轮 B1）后用对应商业化发夹；
   - SNAIL 填入真实 UGI 条码序列替换 22 个 N 占位——**注意**：N 占位期间
     SNAIL 第三级 padlock 特异性比对会静默失效（见
     [丙酰辅酶A基因四方案批量测试报告.md](丙酰辅酶A基因四方案批量测试报告.md)
     的问题 P3），填入真实条码后该检查才生效。

---

## 5. 引擎问题记录（本轮复现，未修改）

本轮在低丰度配置下复现了上一轮报告中的非阻断性问题，**均未修正**：

- **P1 ResourceWarning（160/160 次触发）**：FASTA 句柄未关闭
  （[schemes/common.py `load_first_target`](src/mycoprimer/schemes/common.py#L83-L89)、
  [mining.py `load_fasta`](src/mycoprimer/mining.py#L22-L24) 传文件路径给
  `SeqIO.parse` 未显式关闭）。批量/GUI 长会话的资源泄漏隐患。
- **P3 SNAIL padlock N 占位致第三级特异性检查静默失效**：本轮 SNAIL 仍用
  22 个 N 占位 UGI，padlock 比对恒 0 命中；填入真实条码后恢复。
- **P2 SNAIL 选点间距双重计数**：潜在逻辑问题，本轮无背景过滤、候选更密，
  仍因热力学存活候选稀疏而未表现为产量损失（建议后续修）。
- **P6 Tm 语义不一致**：smFISH 为探针 Tm、HCR3 为 52nt tile Tm、SNAIL 为
  41nt cassette Tm，横向比较需注意。

这些问题不影响本轮"低丰度配置下四方案均可对全部基因产出探针"的结论。

---

## 6. 复现方法

```bash
# 项目目录、Probe 环境
/opt/anaconda3/envs/Probe/bin/python batch_propionyl_lowabundance_test.py
# 汇总：test_data/propionyl_lowabundance_results/propionyl_lowabundance_summary.csv
# 宽表：test_data/propionyl_lowabundance_results/propionyl_lowabundance_wide.csv
```

`src/` 引擎代码在测试前后保持不变（git 状态仅新增测试脚本与结果文件）。
