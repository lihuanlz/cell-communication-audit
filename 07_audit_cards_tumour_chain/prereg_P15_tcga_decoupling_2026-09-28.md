# 预注册卡 P15：致癌 Ras/Raf 突变在病人肿瘤中是否复现"信道解耦"
日期：2026-09-28 状态：FROZEN（取值前）

## 背景
P13 在 MEF 等基因系实测：强致癌 Ras/Raf 突变把 配体→ERK 信道容量压缩 0.3-0.4 bit
（Q61R/BRAF-V600E 最强），且随等位强度分级。该结论目前全部来自培养细胞。
本卡检验其在病人肿瘤中的可观测对应物。

## 逻辑链
若信道被突变卡死（ERK 活性由突变组成性驱动、对上游受体输入失聪），
则在肿瘤组织层面应观察到：pERK 水平与上游受体（EGFR）表达解耦。
对照预言（P12a）：突变动均值不动分布宽度 → 组内 CV 应不变。

## 数据（公开，UCSC Xena / TCGA pan-cancer）
- RPPA 蛋白组：TCGA.RPPA.sampleMap/RPPA_RBN（pERK = MAPK_PT202_Y204；EGFR 蛋白）
- 突变：MC3 nonsilentGene 矩阵（KRAS、BRAF、HRAS、NRAS）
- mRNA：EB++AdjustPANCAN geneExp（EGFR 表达，作蛋白缺失时备份口径）

## 判决线（冻结）
癌种内配对（KRAS 或 BRAF 突变 vs 双野生型），主癌种：
COAD/READ（KRAS/BRAF 高频）、LUAD（KRAS）、SKCM（BRAF/NRAS）、PAAD（KRAS）。
1. 主检验：pERK ~ EGFR 的回归斜率（组内），突变组应显著低于 WT 组
   （交互作用 p < 0.05 且斜率比 < 0.8）。
2. 分级检验：BRAF-V600E（强）vs KRAS-G12X（中）vs 双 WT，斜率应同序递减。
3. 宽度对照：两组 pERK 组内 CV 差异 |ΔCV| < 0.05（P12a 预言的肿瘤版）。
4. 成立 = 主检验过且分级同序；不成立 = 斜率比 ≥ 0.8 或方向反转；不决 = 样本不足
   （任一组 n < 30）。

## 口径登记
- RPPA 为群体组织测量，非单细胞；这是 P13 的"群体投影"，灵敏度低，
  不成立不反杀 P13，成立则强力支持外推。此不对称性取值前登记。
- 混杂：肿瘤纯度、间质污染。用全样本同质过滤，登记纯度分层为次级分析。
