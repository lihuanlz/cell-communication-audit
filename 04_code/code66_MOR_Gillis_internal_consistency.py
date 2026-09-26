# -*- coding: utf-8 -*-
"""
代码66_MOR_Gillis内部一致性审计.py  v1.0.0
============================================================
细胞线4 · GPCR偏向性断裂审计 · 第三受试体：μOR（Gillis et al. 2020, Sci Signal）

目的
----
Gillis 2020 主文表1–3 与 SI 表S1 给出同一批配体×测定单元的四套统计量：
  表1  Emax（% DAMGO）
  表2  pEC50
  表3  log τ（6个部分激动剂 × 5+3 测定）
  表S1 log(τ/KA)（9个配体）
Black–Leff 运算模型在单元格内给出两条独立恒等式（双重冗余）：
  恒等式1  log τ  = log10(Emax/(100−Emax))        （Emax 以系统最大=100 归一）
  恒等式2  log(τ/KA) = log τ + pEC50 − log10(1+τ)
若四表自洽，两条恒等式必须在每个单元格内成立。
另做 检验3：跨测定 pKA = pEC50 − log10(1+τ) 展布（功能亲和力应与测定臂无关——
D2R/代码63 已验证该量在健康数据中跨臂守恒）。

数据：全部由用户手工提供的主文 PDF / SI PDF 经 pdftotext -layout 提取，
逐值核对（主文.txt 行1219–1297；SI.txt 行414–439）。Oxycodone GIRK 为 ND。
注意：SI 表S1 的 GPA 列与 cAMP 列数值与 SEM 完全相同——查原始 PDF 确认
原文即如此（非提取错误），本身即为一处可疑点，如实记录。

立场纪律：本文不做“谁对谁错”的预判。Gillis–Stahl 公案（Stahl 2022 Biochemistry
再分析质疑 Gillis 的 τ 估计）背景下，本审计只交付“不一致地图”。
============================================================
"""
import numpy as np

LIG6 = ["Morphine", "Oxycodone", "Oliceridine", "PZM21", "SR-17018", "Buprenorphine"]
LN10 = np.log(10)

# 每个测定: [v,se]×6 配体, 顺序同 LIG6；None = ND
EMAX = {
 "Nb33":     [71,3, 70,4, 42,8, 38,3, 20,6, 26,3],
 "mGsi":     [80,7, 75,3, 51,7, 52,5, 35,5, 36,3],
 "GPA":      [98,5, 103,4, 85,5, 86,7, 61,13, 79,6],
 "cAMP":     [97,4, 106,7, 86,5, 84,7, 62,12, 86,9],
 "GIRK":     [88,4, None,None, 74,4, 86,3, 78,7, 54,5],
 "GRK2rec":  [66,6, 64,8, 40,3, 36,7, 41,9, 28,5],
 "BarrGRK2": [70,3, 72,5, 58,3, 59,4, 49,7, 34,2]}
PEC50 = {
 "Nb33":     [6.67,.03, 5.93,.12, 7.19,.28, 6.88,.36, 7.48,.46, 8.28,.21],
 "mGsi":     [6.94,.05, 6.22,.11, 7.46,.18, 7.49,.14, 7.15,.33, 8.65,.14],
 "GPA":      [7.72,.16, 6.94,.23, 8.38,.19, 8.16,.19, 7.66,.18, 8.85,.06],
 "cAMP":     [8.10,.26, 7.01,.18, 8.66,.13, 8.64,.25, 7.67,.33, 9.61,.37],
 "GIRK":     [7.98,.10, None,None, 8.43,.12, 8.42,.07, 6.43,.13, 7.97,.21],
 "GRK2rec":  [6.78,.08, 6.28,.09, 7.26,.51, 7.58,.19, 6.99,.19, 7.95,.28],
 "BarrGRK2": [7.31,.13, 6.22,.06, 7.71,.08, 7.56,.07, 6.48,.46, 8.50,.21]}
LOGTAU = {  # 主文表3
 "Nb33":     [0.42,.08, 0.56,.12, -0.34,.12, -0.33,.07, -0.86,.21, -0.62,.09],
 "mGsi":     [0.70,.07, 0.65,.09, -0.12,.05, -0.08,.08, -0.37,.16, -0.40,.05],
 "GPA":      [1.74,.22, 1.86,.26, 1.15,.25, 1.18,.28, 0.66,.37, 0.63,.25],
 "cAMP":     [2.00,.31, 1.87,.22, 1.29,.23, 1.44,.39, 1.04,.28, 1.35,.39],
 "GIRK":     [0.09,.05, None,None, -0.24,.05, -0.18,.04, -0.28,.12, -0.61,.10],
 "GRK2rec":  [0.21,.07, 0.22,.14, -0.30,.18, -0.37,.11, -0.30,.09, -0.57,.08],
 "BarrGRK2": [0.34,.07, 0.35,.09, 0.13,.03, 0.11,.06, -0.04,.15, -0.30,.09]}
LOGR = {    # SI 表S1（cAMP 列与 GPA 列原文完全相同，见文件头注）
 "Nb33":     [6.5,.04, 5.8,.11, 6.7,.17, 6.9,.18, 7.1,.37, 7.6,.2],
 "mGsi":     [6.8,.06, 6.1,.1, 7.4,.08, 7.2,.1, 6.6,.38, 8.0,.14],
 "GPA":      [7.8,.19, 7.0,.24, 8.2,.24, 8.0,.27, 6.1,.58, 8.6,.13],
 "cAMP":     [7.8,.19, 7.0,.24, 8.2,.24, 8.0,.27, 6.1,.58, 8.6,.13],
 "GIRK":     [6.9,.03, None,None, 7.5,.04, 7.3,.03, 5.5,.07, 7.0,.06],
 "BarrGRK2": [7.2,.14, 6.1,.07, 7.5,.07, 7.3,.07, 6.7,.5, None,None]}

def logtau_from_emax(em, eme):
    t = em/(100-em)
    return np.log10(t), (1/LN10)*100*eme/(em*(100-em))

def s_log1ptau(lt, sd_lt):
    t = 10**lt
    return (t/(1+t))*sd_lt   # d log10(1+t)/d log10(t) = t/(1+t)

print("="*100)
print("恒等式1：logτ(表3) vs logτ(由表1 Emax 反演)")
print("="*100)
res1 = []
for a in ["Nb33","mGsi","GPA","cAMP","GIRK","BarrGRK2"]:
    for j,lig in enumerate(LIG6):
        lt, lte = LOGTAU[a][2*j], LOGTAU[a][2*j+1]
        em, eme = EMAX[a][2*j], EMAX[a][2*j+1]
        if lt is None or em is None: continue
        le, sde = logtau_from_emax(em, eme)
        d = lt - le
        z = abs(d)/np.sqrt(lte**2 + sde**2)
        flag = "⚠" if abs(d) > max(0.3, 2*np.sqrt(lte**2+sde**2)) else ""
        res1.append((a, lig, d, z))
        print(f"{a:9s} {lig:14s} 表3={lt:+.2f}±{lte:.2f}  Emax反演={le:+.2f}±{sde:.2f}  Δ={d:+.2f} ({z:.1f}σ) {flag}")
ds = np.array([abs(r[2]) for r in res1]); zs = np.array([r[3] for r in res1])
print(f"\n全体 |Δ| 中位={np.median(ds):.2f} dex, >0.3dex 占比={100*np.mean(ds>0.3):.0f}%, >2σ 占比={100*np.mean(zs>2):.0f}%")

print()
print("="*100)
print("恒等式2：log(τ/KA)(表S1) vs logτ表3 + pEC50 − log10(1+τ表3)")
print("="*100)
res2 = []
for a in ["Nb33","mGsi","GPA","cAMP","GIRK","BarrGRK2"]:
    for j,lig in enumerate(LIG6):
        lt, lte = LOGTAU[a][2*j], LOGTAU[a][2*j+1]
        pe, pee = PEC50[a][2*j], PEC50[a][2*j+1]
        lr, lre = LOGR[a][2*j], LOGR[a][2*j+1]
        if lt is None or lr is None or pe is None: continue
        pred = lt + pe - np.log10(1 + 10**lt)
        sd = np.sqrt(lte**2 + pee**2 + s_log1ptau(lt, lte)**2)
        d = lr - pred
        z = abs(d)/np.sqrt(sd**2 + lre**2)
        flag = "⚠" if abs(d) > max(0.3, 2*np.sqrt(sd**2+lre**2)) else ""
        res2.append((a, lig, d, z))
        print(f"{a:9s} {lig:14s} S1={lr:5.2f}±{lre:.2f}  合成={pred:5.2f}±{sd:.2f}  Δ={d:+.2f} ({z:.1f}σ) {flag}")
ds = np.array([abs(r[2]) for r in res2]); zs = np.array([r[3] for r in res2])
print(f"\n全体 |Δ| 中位={np.median(ds):.2f} dex, >0.3dex 占比={100*np.mean(ds>0.3):.0f}%, >2σ 占比={100*np.mean(zs>2):.0f}%")

print()
print("="*100)
print("检验3：跨测定 pKA = pEC50 − log10(1+τ表3) 展布（应与测定臂无关）")
print("="*100)
A7 = ["Nb33","mGsi","GPA","cAMP","GIRK","GRK2rec","BarrGRK2"]
for j,lig in enumerate(LIG6):
    vals = []
    for a in A7:
        lt, pe = LOGTAU[a][2*j], PEC50[a][2*j]
        if lt is None or pe is None: continue
        lte, pee = LOGTAU[a][2*j+1], PEC50[a][2*j+1]
        pka = pe - np.log10(1 + 10**lt)
        sd = np.sqrt(pee**2 + s_log1ptau(lt, lte)**2)
        vals.append((a, pka, sd))
    arr = np.array([v[1] for v in vals])
    print(f"\n{lig}: 展布={arr.max()-arr.min():.2f} dex")
    for a, pka, sd in vals:
        print(f"   {a:9s} pKA={pka:5.2f}±{sd:.2f}")
    best = (0, None)
    for i in range(len(vals)):
        for k in range(i+1, len(vals)):
            z = abs(vals[i][1]-vals[k][1])/np.sqrt(vals[i][2]**2+vals[k][2]**2)
            if z > best[0]:
                best = (z, (vals[i][0], vals[k][0], vals[i][1]-vals[k][1]))
    print(f"   最大两两不一致: {best[0]:.1f}σ ({best[1][0]} vs {best[1][1]}, Δ={best[1][2]:+.2f})")
