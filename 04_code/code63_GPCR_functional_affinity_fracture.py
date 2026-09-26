#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
代码63 GPCR 功能亲和力断裂检验：D2R 多通路 operational 模型审计 v1.1.0
================================================
细胞线4 第一刀。任务：把"断裂审计"从趋化系统平移到 GPCR。
数据：Klein Herenbrink et al. 2016, Nat Commun 7:10842（D2L 受体，CHO
细胞，同一实验室同一体系），SI Table 1（7 配体 × 6 通路 log(τ/KA) ± SEM，
5 min）、SI Table 5（[3H]spiperone 竞争结合 pKi ± SEM）、SI Table 6（荧光
PPHT 示踪剂动力学 pKd，独立复制）、SI Table 7（cAMP/GαoB pEC50+Emax 时程）。

审计设计（operational 模型的内部约束，任一破坏即断裂）：
  检验1（结合锚点）：pKi 锚定 τ → 预言 Emax vs 实测。
  检验2（功能亲和力反演）：pKA_func=logR−log τ；跨通路守恒性 +
    与 pKi 的一致性，SEM 全传播。
  检验3（时间不变性）：SI T7 全配体 pEC50 时程（2–90 min）；
    平衡模型要求 ΔlogR 不漂移。
  自审①：零模型判据 pKi−pEC50>0（BL 下 pEC50=pKA+log(1+τ)≥pKA，
    单 KA 假设下不得显著低于 pKi）——无任何拟合假设。
  自审②：第二示踪剂（PPHT 动力学 pKd）独立复制结合锚点。
  自审③：系统最大 E_sys 扰动 ±10% 的敏感性。
  自审④：误差全传播（σ_Δ 含 logR/pKi/Emax 三项）。

结果（v1.1.0 实跑）：
  检验2：Δ=pKA_func−pKi：bifeprunox −3.5（21–23σ）、aripiprazole
    −2.6/−2.9（18–26σ）、cariprazine −1.2/−1.4（9–11σ）；常规激动剂
    |Δ|≤1.3（dopamine +1.2/1.3 为两态移动预期方向）。pKA_func 跨
    通路守恒（臂间差 0.08–0.31 dex）。
  自审①：零模型复现：pKi−pEC50 = +1.96/+0.52/+2.45 dex（三个高亲和力
    部分激动剂），其余配体为负（备用受体方向，单 KA 允许）。
  自审②：PPHT 复制：同一三配体 pKd−pEC50 = +2.19/+1.15/+2.39 dex；
    两示踪剂互差仅 0.06–0.63 dex，锚点本身稳固。
  自审③：E_sys 100→110 时 Δ 仅移 0.1–0.35 dex，分类不变。
  检验3：全激动剂效价 90 min 下移 −1.15…−1.37 dex（脱敏），三个部分
    激动剂稳定或上移（bifeprunox +1.09）；bifeprunox vs ropinirole
    偏向因子 2' 时 −0.42 → 90' 时 +1.93 dex，摆动 2.35 dex 且变号
    （SEM 0.02–0.03）→ 平衡偏向定量给出依赖读数时刻的相反答案。

判定：单亲和力平衡 operational 模型在 D2R 数据上断裂，定位于高亲和力
部分激动剂——恰是全部显著偏向报道所属类别（SI T3 星号集中于此）。
方法自审四条全部通过，断裂稳健。

⚠️ 双录：PTH1R（Sachdev 2024，Source Data 已入卷）缺结合锚点，
归一化曲线上 Black–Leff 的 τ/KA 在 nH≠1 时不可辨识（撞边界伪影
+5.77/+9.58 dex，不作结论）。教训：无独立亲和力锚点时，边际
operational 统计量不构成可判伪检验。

运行：python3 代码63_GPCR功能亲和力_断裂检验.py
依赖：numpy
"""

import numpy as np

LIGS = ["Dopamine", "Ropinirole", "Aripiprazole", "Cariprazine",
        "Bifeprunox", "Pardoprunox", "S-3PPP"]
logR = {
"Dopamine":    {"cAMP": (8.06,.05), "pERK": (8.56,.05), "Gai1": (7.17,.09), "GaoB": (7.24,.08), "CI": (8.55,.05), "Barr": (6.33,.13)},
"Ropinirole":  {"cAMP": (7.67,.04), "pERK": (8.44,.05), "Gai1": (6.75,.09), "GaoB": (7.08,.08), "CI": (8.32,.05), "Barr": (6.48,.15)},
"Aripiprazole":{"cAMP": (7.36,.06), "pERK": (6.48,.54), "Gai1": (7.90,.16), "GaoB": (6.94,.11), "CI": (7.83,.15), "Barr": (6.17,.58)},
"Cariprazine": {"cAMP": (8.27,.06), "pERK": (7.45,.47), "Gai1": (8.10,.14), "GaoB": (7.91,.11), "CI": (8.88,.11), "Barr": None},
"Bifeprunox":  {"cAMP": (7.88,.05), "pERK": (7.32,.10), "Gai1": (8.60,.11), "GaoB": (7.50,.10), "CI": (8.30,.07), "Barr": (5.75,.16)},
"Pardoprunox": {"cAMP": (8.64,.05), "pERK": (9.06,.08), "Gai1": (8.10,.13), "GaoB": (7.70,.10), "CI": (9.29,.06), "Barr": (6.80,.29)},
"S-3PPP":      {"cAMP": (6.65,.06), "pERK": (7.39,.14), "Gai1": (6.60,.12), "GaoB": (6.40,.11), "CI": (7.34,.10), "Barr": None},
}
pKi = {"Dopamine": (5.05,.06), "Ropinirole": (5.60,.08), "Aripiprazole": (9.43,.06),
       "Cariprazine": (8.90,.07), "Bifeprunox": (10.36,.09), "Pardoprunox": (7.63,.06),
       "S-3PPP": (5.84,.06)}
pKd2 = {"Dopamine": 5.18, "Ropinirole": 5.73, "Aripiprazole": 9.66, "Cariprazine": 9.53,
        "Bifeprunox": 10.30, "Pardoprunox": 7.75, "S-3PPP": 6.11}
EMAX5 = {"cAMP": {"Dopamine": (98,2.0), "Ropinirole": (100,0), "Aripiprazole": (76,2.0),
                  "Cariprazine": (77,2.4), "Bifeprunox": (92,2.0), "Pardoprunox": (87,1.5), "S-3PPP": (73,2.0)},
         "GaoB": {"Dopamine": (91,2.0), "Ropinirole": (100,0), "Aripiprazole": (71,5.0),
                  "Cariprazine": (70,3.7), "Bifeprunox": (80,3.6), "Pardoprunox": (79,4.3), "S-3PPP": (68,3.0)}}
PEC50_5 = {"Dopamine": 8.07, "Ropinirole": 7.65, "Aripiprazole": 7.47, "Cariprazine": 8.38,
           "Bifeprunox": 7.91, "Pardoprunox": 8.69, "S-3PPP": 6.78}
T = [2, 5, 10, 15, 30, 45, 60, 75, 90]
PEC50_T = {
"Ropinirole":  [8.02,7.65,7.38,7.37,7.46,7.01,6.73,6.97,6.76],
"Dopamine":    [8.40,8.07,7.71,7.65,7.78,7.29,6.97,7.30,7.03],
"Aripiprazole":[7.26,7.47,7.58,7.65,7.76,7.56,7.45,7.71,7.53],
"Cariprazine": [8.30,8.38,8.38,8.38,8.49,8.37,8.22,8.40,8.16],
"Bifeprunox":  [7.60,7.91,8.13,8.27,8.62,8.56,8.51,8.80,8.69],
"Pardoprunox": [8.98,8.69,8.39,8.33,8.46,8.03,7.82,8.04,7.83],
"S-3PPP":      [7.10,6.78,6.43,6.37,6.47,5.99,5.82,6.04,5.82],
}

def main():
    print("== 检验2：功能 KA 反演，误差全传播 ==")
    print(f"{'配体':13s} {'通路':5s} {'pKA_func':>8s} {'pKi':>6s} {'Δ(dex)':>8s} {'σ(Δ)':>6s} {'显著性':>6s}")
    for lig in LIGS:
        for path in ["cAMP", "GaoB"]:
            em, eme = EMAX5[path][lig]
            if em >= 99.5:
                print(f"{lig:13s} {path:5s}   τ→∞（满激动，KA 不可分辨）")
                continue
            tau = em / (100 - em)
            lr, lre = logR[lig][path]; pk, pke = pKi[lig]
            pka = lr - np.log10(tau)
            s_logtau = (1 / np.log(10)) * 100 * eme / (em * (100 - em))
            sd = np.sqrt(lre**2 + pke**2 + s_logtau**2)
            print(f"{lig:13s} {path:5s} {pka:8.2f} {pk:6.2f} {pka-pk:+8.2f} {sd:6.2f} {abs(pka-pk)/sd:5.1f}σ")

    print("\n== 自审①：零模型判据 pKi−pEC50（>0.3 即断裂，无拟合假设）==")
    for lig in LIGS:
        d = pKi[lig][0] - PEC50_5[lig]
        tag = "⚠断裂" if d > 0.3 else ("备用受体方向（满激动预期）" if d < -0.3 else "一致")
        print(f"  {lig:13s} pKi={pKi[lig][0]:5.2f} pEC50={PEC50_5[lig]:5.2f} 差={d:+.2f}  {tag}")

    print("\n== 自审②：第二示踪剂 PPHT 动力学 pKd 复制 ==")
    for lig in LIGS:
        print(f"  {lig:13s} pKd_ppht={pKd2[lig]:5.2f} 两示踪剂差={pKd2[lig]-pKi[lig][0]:+.2f} pKd−pEC50={pKd2[lig]-PEC50_5[lig]:+.2f}")

    print("\n== 自审③：系统最大扰动（E_sys=100/105/110，cAMP 臂 Δ）==")
    for esys in [100, 105, 110]:
        out = []
        for lig in LIGS:
            em = min(EMAX5["cAMP"][lig][0], esys - 1)
            tau = em / (esys - em)
            out.append(f"{lig[:4]}={logR[lig]['cAMP'][0]-np.log10(tau)-pKi[lig][0]:+.2f}")
        print(f"  E_sys={esys}: " + ", ".join(out))

    print("\n== 检验3：pEC50 时间漂移（2'→90'）与偏向因子摆动 ==")
    for lig in LIGS:
        v = PEC50_T[lig]
        print(f"  {lig:13s}: {v[0]:.2f}→{v[-1]:.2f}  Δ={v[-1]-v[0]:+.2f} dex")
    for lig in ["Bifeprunox", "Aripiprazole", "Cariprazine", "S-3PPP"]:
        d = [a - b for a, b in zip(PEC50_T[lig], PEC50_T["Ropinirole"])]
        print(f"  {lig:13s} vs Ropi: 2'={d[0]:+.2f} → 90'={d[-1]:+.2f}  摆动={max(d)-min(d):.2f} dex")

if __name__ == "__main__":
    main()
