# 判词卡：Track C2 —— 真实 (α, ν) 签名落判别图

**日期：2026-09-18 ｜ 性质：定理→数据重分析（拓扑选择定理的实测验证）**
**谱系：p53 备忘录 §7 步骤 4 / TCS项目分析报告 Track C2（开放两年）→ 本次关闭（文本级+数字化级，非逐像素抠图级）。**

## 0. 一句话判词

**真实数据的 (α, ν) 签名与拓扑选择定理的判别区完全咬合：DSB 双臂深落 EXC 区（α≈0, ν>0），UV 落模拟区（α≈1, ν≈0）——同一蛋白同一细胞系（MCF7）内签名-拓扑共变，定理从"数值验证"升级为"实测验证"。**

## 1. 数字（全部带误差棒与来源层级）

| 数据 | α = ∂lnA/∂lnD | ν = ∂N/∂lnD | 落区 | 来源 |
|---|---|---|---|---|
| NCS（DSB） | **−0.001 ± 0.113**（95%CI [−0.23, 0.22]） | 2.0 ± 0.5 | **EXC** | α：仓内数字化 Fig1G bootstrap×20000；ν：Mönke 2017 经 TCS 报告 |
| γ-IR（DSB） | 0.05 ± 0.10 | 1.4 ± 0.4 | **EXC** | Batchelor 2011 原文"幅度与剂量无关"（仓内 fulltext.xml 锚）+ Lahav 2004 N:1→~6/0.1–10 Gy |
| UV | **0.981 ± 0.114**（95%CI [0.78, 1.22]） | 0.0 +0.2 | **模拟/NF 签名** | α：仓内数字化 Fig1H bootstrap×20000；ν：原文"单脉冲" |
| 模型锚 NF | 1.13（备忘录模型 A 表） | ≈0 | — | 08-02 备忘录存档值 |
| 模型锚 EXC | 0.04 | τ_r/T（参数） | — | 代码51/备忘录 |

**两区分离度**：α_NCS 与 α_UV 相距 ≈ 8.6 个联合标准差——判别对测量噪声鲁棒，定理第三条（数量级分离）实测成立。

## 2. 原文锚（仓内 fulltext.xml，直接引用）

- "DSBs trigger a series of p53 pulses with **fixed amplitude and duration, independent of the damage dose**, whereas **the number of pulses increases with higher damage**."
- UV 臂："single pulse that increases in amplitude and duration **in proportion to the UV dose**."
- GZ2006（本次抓取全文复核）："mean amplitude and period … **did not appear to significantly depend on irradiation level**"（Suppl. Fig S5）；振幅 CV≈70% vs 周期 <20%。

## 3. 副产品：计数律斜率反演修复时标

引理 3（N ≈ (τ_r/T)·ln(D0/D_c)）的斜率即修复时间常数的可辨识估计：
**τ_r ≈ ν·T = 1.4 × 5.5 h ≈ 7.7 h**——与 DSB 修复快分量时标（小时量级）一致。计数律斜率从"拟合参数"升格为"测量量"。

## 4. 诚实边界（不越界声明）

1. 本次为**文本级 + 已数字化级**重分析：α 值来自我们自己的数字化（Fig1G/H，bootstrap 误差棒），ν 值来自论文原文陈述与文献文本值（Lahav N:1→6、Mönke ν≈2 估计）——**未做逐像素抠图**。若需 ν 的独立误差棒（而非引用宽区间），需 WebPlotDigitizer 抠 Lahav 2004 Fig/Mönke 2017 Fig1D/E，留作升级项；
2. NCS 的 ν=2.0±0.5 引用自 TCS项目分析报告对 Mönke 2017 的估计，误差棒是保守区间不是拟合产物；
3. γ-IR 的 α=0.05±0.10 中 0.05 来自对账卡对文献 CV 的折算，±0.10 为保守区间；核心结论（α≈0）有原文直接陈述背书；
4. 判别图模型锚 NF 用 08-02 备忘录存档表值，未重跑（参数套以备忘录为准）。

## 5. 产物

- 代码：`03_细胞线3/代码与图/代码52_TrackC2_判别图.py`
- 图：`03_细胞线3/结果/代码52_TrackC2_alpha-nu判别图.png`
- 数据：`03_细胞线3/结果/代码52_TrackC2_结果.json`

## 6. 结论

**Track C2 关闭（文本级+数字化级）。** 升级项（逐像素抠图出 ν 独立误差棒）不影响任何现有结论，仅收紧 ν 的不确定性。
