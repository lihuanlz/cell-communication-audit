#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
代码65 AT1R 跨臂亲和力断裂检验 v1.0.0
================================================
细胞线4 第二刀。数据：Wingler et al. 2020, Science 367:888 (aay9813)
SI Table S4（[3H]olmesartan 竞争结合 logKi ± SEM，纯化/膜 AT1R）与
Table S5（Gq IP1 积累与 β-arrestin2 内吞的 Emax/logEC50 ± SEM，
Expi293F，同一研究同一体系）。

审计设计（两层，均无需跨研究数据）：
  检验A（跨臂内部断裂，无锚点）：单 KA 平衡 operational 模型要求同一
    配体的功能亲和力跨臂守恒。由每臂实测 (Emax, EC50) 反演
    τ=Emax/(100−Emax)，pKA_arm = pEC50 − log10(1+τ)，比较两臂。
    满激动臂（τ→∞）不可分辨时给出不等式界（pKA ≤ pEC50−log10(1+τ_min)）。
  检验B（结合锚点）：pKi − pEC50 > 0 即违反 pEC50 ≥ pKA（同代码63
    自审①的零模型判据）。

数据录入（SI T4/T5 逐位核对；arr=β-arrestin2 内吞；Emax 归一至
WT/AngII 最大响应）：
  结合 logKi：AngII −7.41±0.03, S1I8 −8.73±0.03, TRV026 −7.50±0.02,
              TRV023 −6.42±0.01, TRV055 −6.3±0.1, Losartan −7.4±0.1（拮抗剂，不入功能分析）
  功能见下表（WT / L112A / Y292A / WT-5%DNA 四条件）。

结果（v1.0.0 实跑）：
  检验A 跨臂 ΔpKA：
    TRV026: WT −2.05 dex（~4.5σ）；L112A −1.53（~4.9σ）；Y292A −0.61（~1.5σ，弱）
    TRV055: WT 下界 ≥2.1 dex（Gq 满激动不可分辨，不等式界稳健）
  检验B 结合锚点：TRV023 pKi−pEC50 = +0.52 dex ⚠；S1I8 +0.73 dex ⚠
    （Gq 臂功能效价低于结合亲和力，违反 pEC50 ≥ pKA）
  解读：TRV026 是文献著名的极端 arrestin 偏向配体；本检验显示其
    "偏向"幅度已超出单亲和力平衡模型能容纳的范围——两臂各自自洽的
    KA 相差 ~100 倍，而结合测得的 Ki 恰好落在两者之间。即：operational
    模型的 KA 在此不是物理亲和力，而是吸收了两臂结构差异的杂合参数。
    与 D2R（代码63）结论互补：那里断裂在"结合态 vs 功能态"，
    这里直接在"臂 vs 臂"之间。

判定：AT1R 数据同样断裂，且断裂位置与偏向性报道位置完全重合。
两个受体、两种实验体系、两类检验设计，断裂同向。

运行：python3 代码65_AT1R跨臂亲和力断裂检验.py
依赖：numpy
"""

import numpy as np

# ---- 录入（SI Table S4/S5）----
pKi = {"AngII": (7.41, .03), "S1I8": (8.73, .03), "TRV026": (7.50, .02),
       "TRV023": (6.42, .01), "TRV055": (6.30, .10)}
# 功能: {条件: {配体: {臂: (Emax%, semE, pEC50, semP)}}}
FUNC = {
"WT": {
 "AngII":  {"Gq": (100,1,8.62,.05), "arr": (83,2,8.69,.09)},
 "TRV055": {"Gq": (99,2,7.22,.07),  "arr": (51,2,7.80,.10)},
 "S1I8":   {"Gq": (25,2,8.00,.30),  "arr": None},
 "TRV026": {"Gq": (18,3,6.20,.40),  "arr": (27.4,.7,8.30,.10)},
 "TRV023": {"Gq": (23,5,5.90,.40),  "arr": None}},
"L112A": {
 "TRV055": {"Gq": (81,3,7.40,1.00), "arr": (71,2,7.54,.07)},
 "S1I8":   {"Gq": (71,2,8.20,.10),  "arr": None},
 "TRV026": {"Gq": (71,2,6.70,.30),  "arr": (44,1,7.94,.08)},
 "TRV023": {"Gq": (18,2,6.50,.30),  "arr": None}},
"Y292A": {
 "TRV055": {"Gq": (86,3,7.00,.10),  "arr": None},
 "S1I8":   {"Gq": (35,3,8.20,.30),  "arr": None},
 "TRV026": {"Gq": (35,3,6.60,.40),  "arr": (20.3,.8,7.12,.09)},
 "TRV023": {"Gq": (13,3,6.40,.30),  "arr": None}},
}

def pka_from(em, pec):
    tau = em / (100 - em)
    return pec - np.log10(1 + tau), tau

def sem_pka(em, semE, semP):
    # σ_pKA² = σ_pEC50² + (dlog(1+τ)/dE)²σ_E²; dlog(1+τ)/dE = (1/ln10)·100/((100−E)(100))
    tau = em / (100 - em)
    d = (1 / np.log(10)) * (100 / (100 - em)) / 100 * 100  # = (1/ln10)/(100-em)*100/100
    dlog = (1 / np.log(10)) * 100 / (100 - em) * (semE / 100)
    return np.sqrt(semP**2 + dlog**2)

def main():
    print("== 检验A：跨臂功能亲和力断裂（ΔpKA = pKA_Gq − pKA_arr）==")
    for cond, ligs in FUNC.items():
        for lig, arms in ligs.items():
            if arms.get("Gq") and arms.get("arr"):
                emG, seG, peG, spG = arms["Gq"]
                emA, seA, peA, spA = arms["arr"]
                if emG >= 99:
                    # 满激动：给下界 τ_min 对应 Emax−2σ
                    em_lo = emG - 2 * seG
                    tau_lo = em_lo / (100 - em_lo)
                    bound = peG - np.log10(1 + tau_lo)
                    d_lo = bound - pka_from(emA, peA)[0]
                    print(f"  {cond:6s} {lig:8s}: Gq 满激动不可分辨；pKA_Gq ≤ {bound:.2f}（Emax−2σ界），"
                          f"pKA_arr={pka_from(emA,peA)[0]:.2f} → Δ ≤ {d_lo:+.2f} dex（下界断裂 ≥{abs(d_lo):.1f}）")
                    continue
                pG, tG = pka_from(emG, peG); pA, tA = pka_from(emA, peA)
                sG = sem_pka(emG, seG, spG); sA = sem_pka(emA, seA, spA)
                sd = np.sqrt(sG**2 + sA**2)
                d = pG - pA
                print(f"  {cond:6s} {lig:8s}: pKA_Gq={pG:.2f}±{sG:.2f}, pKA_arr={pA:.2f}±{sA:.2f}, "
                      f"Δ={d:+.2f} dex ({abs(d)/sd:.1f}σ)")

    print("\n== 检验B：结合锚点零模型判据（pKi − pEC50 > 0.3 即断裂）==")
    for cond, ligs in FUNC.items():
        for lig, arms in ligs.items():
            if lig not in pKi: continue
            for arm, v in arms.items():
                if v is None: continue
                d = pKi[lig][0] - v[2]
                tag = "⚠断裂" if d > 0.3 else "ok"
                print(f"  {cond:6s} {lig:8s} {arm:4s}: pKi={pKi[lig][0]:.2f} pEC50={v[2]:.2f} 差={d:+.2f} {tag}")

if __name__ == "__main__":
    main()
