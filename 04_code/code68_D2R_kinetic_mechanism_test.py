# -*- coding: utf-8 -*-
"""
代码68_D2R动力学机制检验.py  v1.0.0
============================================================
细胞线4 · 挖深轮：D2R 断裂的动力学机制裁决

问题：代码63 判定 bifeprunox/aripiprazole/cariprazine 三个高亲和力部分激动剂
存在结合态(pKi≈8.9–10.4)与功能态(pKA_func≈6.9–7.8)的 −1.2…−3.5 dex 断裂，
且偏向因子随时间反转（bifeprunox vs ropinirole −0.42 → +1.93 dex）。
Klein Herenbrink 2016 自己提出"动力学语境"定性解释（慢解离配体未达平衡）。
本代码做定量裁决：用该文 SI 表5/表6 实测的 kon/koff，构造非平衡占用率
  ρ(L,t) = [kon·L/(kon·L+koff)]·(1−e^{−(kon·L+koff)·t})
代入运算模型（n=1，E = Em·τρ/(1+ρ(τ−1))），τ 由 90' Emax 反演（Esys=100，
ropinirole 满激动剂参照），生成表观浓度-反应曲线并提取表观 pEC50(t)，
与 SI 表7 实测的 8 个时间点 pEC50 轨迹对账。**零自由参数**（Kd 不拟合，
直接取动力学测量值）。

判决逻辑：非平衡只能把表观 potency 暂时压到平衡值以下；t→∞ 时模型地板为
pKd+log10(1+τ) ≥ pKd。若观测 pEC50 在任何时刻都低于该地板 1 dex 以上，
动力学解释在定量上不成立，断裂存活。

数据（Klein Herenbrink 2016, Nat Commun 7:10842, SI）：
  表5：[3H]spiperone 示踪竞争结合动力学（仅3个慢配体可测）
  表6：PPHT-red 荧光示踪 Tag-lite 结合动力学（7配体）
  表7：cAMP/Gαo/CI × 2–90 min × 7配体 的 pEC50/Emax
  ——注意：表5 与表6 对 cariprazine 的 kon 差 82 倍（pKd 7.50 vs 9.56），
    论文自身结合数据内部不一致，如实记录并做双示踪剂敏感性分析。
============================================================
"""
import numpy as np
from scipy.interpolate import interp1d

LIGS7 = ["Ropinirole","Dopamine","Aripiprazole","Cariprazine","Bifeprunox","Pardoprunox","S-3PPP"]
# SI 表6（PPHT）
KON6 = {"Ropinirole":1.46e6,"Dopamine":3.14e5,"Aripiprazole":1.01e9,"Cariprazine":1.27e9,
        "Bifeprunox":1.84e8,"Pardoprunox":1.25e8,"S-3PPP":3.25e6}
KOFF6= {"Ropinirole":2.60,"Dopamine":2.00,"Aripiprazole":0.21,"Cariprazine":0.35,
        "Bifeprunox":0.01,"Pardoprunox":2.28,"S-3PPP":1.51}
# SI 表5（spiperone，仅三慢配体）
KON5 = {"Aripiprazole":1.31e8,"Cariprazine":1.55e7,"Bifeprunox":1.07e8}
KOFF5= {"Aripiprazole":0.14,"Cariprazine":0.49,"Bifeprunox":0.01}
PKI  = {"Aripiprazole":9.43,"Cariprazine":8.90,"Bifeprunox":10.36}  # 表5 pKi

TIMES=[2,5,10,15,30,45,60,75,90]
T7={
"cAMP":{
 "Ropinirole":([8.02,7.65,7.38,7.37,7.46,7.01,6.73,6.97,6.76],[100]*9),
 "Dopamine":([8.40,8.07,7.71,7.65,7.78,7.29,6.97,7.30,7.03],[98,98,98,99,102,101,105,98,102]),
 "Aripiprazole":([7.26,7.47,7.58,7.65,7.76,7.56,7.45,7.71,7.53],[75,76,75,72,71,64,54,71,65]),
 "Cariprazine":([8.30,8.38,8.38,8.38,8.49,8.37,8.22,8.40,8.16],[77,77,75,72,72,66,55,72,67]),
 "Bifeprunox":([7.60,7.91,8.13,8.27,8.62,8.56,8.51,8.80,8.69],[89,92,91,89,89,88,81,90,88]),
 "Pardoprunox":([8.98,8.69,8.39,8.33,8.46,8.03,7.82,8.04,7.83],[87,87,86,87,85,82,74,84,81]),
 "S-3PPP":([7.10,6.78,6.43,6.37,6.47,5.99,5.82,6.04,5.82],[73,73,70,70,71,64,53,66,60])},
"Gao":{
 "Ropinirole":([6.97,7.02,7.17,7.33,7.51,7.38,7.29,7.40,7.36],[100]*9),
 "Dopamine":([7.26,7.31,7.40,7.49,7.30,7.35,7.35,7.47,7.33],[93,91,97,95,101,100,101,99,95]),
 "Aripiprazole":([6.91,7.07,7.36,7.60,7.64,7.72,7.68,7.39,7.49],[80,71,62,63,68,60,60,87,69]),
 "Cariprazine":([7.97,8.05,8.14,8.25,8.36,8.55,8.58,8.71,8.71],[71,70,67,65,67,61,60,52,63]),
 "Bifeprunox":([7.28,7.57,7.89,8.12,8.39,8.52,8.49,8.66,8.67],[94,80,77,79,87,88,80,94,90]),
 "Pardoprunox":([7.78,7.77,8.00,8.21,8.70,8.71,8.64,8.34,8.58],[82,79,76,81,84,85,86,80,87]),
 "S-3PPP":([6.46,6.54,6.57,6.67,7.02,6.78,6.79,6.62,6.77],[64,68,73,76,74,72,72,68,77])},
"CI":{
 "Ropinirole":([8.33,8.29,8.10,7.97,7.77,7.71,7.62,7.55,7.46],[100]*9),
 "Dopamine":([8.55,8.58,8.36,8.25,8.05,8.03,7.97,7.91,7.86],[95,96,101,103,105,106,106,106,106]),
 "Aripiprazole":([7.65,8.34,8.62,8.80,8.74,9.09,9.17,9.28,9.37],[27,30,25,24,22,21,19,18,17]),
 "Cariprazine":([8.94,9.28,9.38,9.45,9.50,9.52,9.55,9.57,9.60],[33,39,34,35,33,33,32,31,30]),
 "Bifeprunox":([7.79,8.46,8.80,9.04,9.26,9.47,9.58,9.64,9.70],[61,67,57,58,53,51,49,46,44]),
 "Pardoprunox":([9.30,9.41,9.22,9.20,9.06,9.02,8.95,8.87,8.80],[75,76,70,67,62,59,57,54,52]),
 "S-3PPP":([7.66,7.67,7.58,7.58,7.50,7.41,7.39,7.40,7.41],[46,46,36,32,26,23,21,19,17])}}

Lgrid = np.logspace(-11, -4.5, 400)
xgrid = np.log10(Lgrid)

def apparent(lig, t, tau, Em, kon, koff):
    kobs = kon*Lgrid + koff
    req  = Lgrid*kon/kobs
    r    = req*(1-np.exp(-kobs*t))
    y    = Em*tau*r/(1+r*(tau-1))
    em   = y[-1]; half = 0.5*em
    if not np.all(np.isfinite(y)) or y.max() < half or y[0] > half:
        return np.nan
    return -float(np.interp(half, y, xgrid))   # pEC50 = -log10(EC50)

def run(KON, KOFF, tag, ligands):
    print(f"\n########## 示踪剂 {tag} ##########")
    for a in ["cAMP","Gao","CI"]:
        print(f"--- {a} ---  (残差 = model − obs, dex)")
        for lig in ligands:
            em90 = T7[a][lig][1][8]
            tau  = em90/(100-em90) if em90 < 99.5 else 200.0
            pred = [apparent(lig, t, tau, 100, KON[lig], KOFF[lig]) for t in TIMES]
            obs  = T7[a][lig][0]
            pKd  = -np.log10(KOFF[lig]/KON[lig])
            dm   = pred[8]-pred[0] if np.isfinite(pred[0]) and np.isfinite(pred[8]) else np.nan
            print(f"{lig:14s} pKd={pKd:5.2f} τ={tau:6.1f}  漂移 model={dm:+.2f}/obs={obs[8]-obs[0]:+.2f}"
                  f"  90'残差={(pred[8]-obs[8]):+.2f}")
            print("    残差 2'→90': " + " ".join(f"{(p-o):+.2f}" if np.isfinite(p) else "  nan"
                                                 for p,o in zip(pred,obs)))

print("表5 vs 表6 动力学亲和自对比：",
      {l: f"pKd5={-np.log10(KOFF5[l]/KON5[l]):.2f} / pKd6={-np.log10(KOFF6[l]/KON6[l]):.2f} / pKi={PKI[l]}"
       for l in KON5})
run(KON6, KOFF6, "表6（PPHT-red）", LIGS7)
run(KON5, KOFF5, "表5（[3H]spiperone）", list(KON5))
