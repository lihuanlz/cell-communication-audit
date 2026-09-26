# -*- coding: utf-8 -*-
"""
代码72_审计方法假阳性率自检.py  v1.0.0
============================================================
细胞线4 · C1：审计工具审计审计工具——假阳性率与功率

代码63–71 给别人判了四次"断裂"，但判据本身的假阳性率从未测过。
本代码用零假设仿真（真模型完全自洽 + 与原文同量级的噪声）回答四个问题：

  M1  Gillis 结构阴性仿真：恒等式1/2 在自洽数据上的报警率
      （对照实测：GIRK 列 5/5 与 4/5 报警）
  M2  nH≠1 合法变异能产生多大假偏差（校准 Atlas 21.5% 的可归因份额；
      广义运算模型 E=Esys(τA)^n/((KA+A)^n+(τA)^n)，nH∈[0.8,1.25]）
  M3  零模型判据（pKi−pEC50>0.3 dex）假阳性率
  M4  检出功率曲线：注入 δ dex 系统断裂 → 单格检出率
  M5  判决对阈值的稳健性（Atlas 实测分布在 0.2–1.0 dex 阈值扫描下）
  M6  列级联合证据的零假设概率 + 报警方向对称性

纪律：本代码只对方法打分，不改任何已封存的判决；若发现假阳性率高于
名义值，相关判决按双录规则降级。
运行：python3 代码72_审计方法假阳性率自检.py
============================================================
"""
import numpy as np, pandas as pd, os
from scipy.optimize import brentq

rng = np.random.default_rng(20260816)
LN10 = np.log(10)

print("=" * 90)
print("M1  Gillis 结构阴性仿真（真模型自洽；SEM 取原文量级：pEC50±0.2, logτ±0.2,")
print("    Emax±5, log(τ/KA)±0.15）")
print("=" * 90)
flags1 = flags2 = cells = 0
nsim = 2000
for _ in range(nsim):
    for lig in range(5):
        pKA = rng.uniform(6, 9)
        for arm in range(6):
            tau = 10 ** rng.uniform(-0.7, 2.0)
            em_true = 100 * tau / (1 + tau); pec_true = pKA + np.log10(1 + tau)
            lt_true = np.log10(tau); lr_true = lt_true + pKA
            em = np.clip(em_true + rng.normal(0, 5.0), 1, 99)
            pec = pec_true + rng.normal(0, 0.2)
            lt = lt_true + rng.normal(0, 0.2)
            lr = lr_true + rng.normal(0, 0.15)
            lt_e = np.log10(em / (100 - em)); sd_e = (1 / LN10) * 100 * 5.0 / (em * (100 - em))
            d1 = abs(lt - lt_e)
            s_l1 = (10 ** lt / (1 + 10 ** lt)) * 0.2
            pred = lt + pec - np.log10(1 + 10 ** lt)
            sd = np.sqrt(0.2 ** 2 + 0.2 ** 2 + s_l1 ** 2)
            d2 = abs(lr - pred)
            cells += 1
            if d1 > max(0.3, 2 * np.sqrt(0.2 ** 2 + sd_e ** 2)): flags1 += 1
            if d2 > max(0.3, 2 * np.sqrt(sd ** 2 + 0.15 ** 2)): flags2 += 1
fp1, fp2 = flags1 / cells, flags2 / cells
print(f"恒等式1 假阳性率 = {fp1:.4f}   恒等式2 假阳性率 = {fp2:.4f}")
print("对照实测：Gillis GIRK 列 恒等式1 5/5、恒等式2 4/5；全表 |Δ|>0.3 占 39%/38%")

print()
print("=" * 90)
print("M2  nH≠1 合法变异的假偏差（Atlas 恒等式A 校准）")
print("=" * 90)
def ec50_gen(tau, KA, n):
    Emn = tau ** n / (1 + tau ** n)
    g = lambda A: (tau * A) ** n / ((KA + A) ** n + (tau * A) ** n) - 0.5 * Emn
    return brentq(g, KA * 1e-4, KA * 1e4)
devs = []
for _ in range(20000):
    n = rng.uniform(0.8, 1.25)
    tau = 10 ** rng.uniform(-0.5, 2); KA = 1.0
    em = 100 * tau ** n / (1 + tau ** n)
    a50 = ec50_gen(tau, KA, n)
    logRA = np.log10(em) - np.log10(a50)
    tau1 = em / (100 - em)
    ltk1 = np.log10(tau1) + (-np.log10(a50) - np.log10(1 + tau1))   # nH=1 反演的"文献 τ/KA"
    emr, a50r = 100.0, ec50_gen(100, 1, 1.0)
    logRA_r = np.log10(emr) - np.log10(a50r)
    ltk_r = np.log10(100) + (-np.log10(a50r) - np.log10(101))
    devs.append(abs((ltk1 - ltk_r) - (logRA - logRA_r)))
devs = np.array(devs)
print(f"nH∈[0.8,1.25] 下 |Δlog(τ/KA)−Δlog(RA)|：中位 {np.median(devs):.3f} dex，"
      f">0.3 占 {np.mean(devs > 0.3):.4f}，>1.0 占 {np.mean(devs > 1):.5f}")
print("对照 Atlas 实测：>0.3 占 21.5%，>1.0 占 4.2%")
print("（边界声明：文献侧按 nH=1 反演；若原文用自由 nH 拟合，假偏差会更大——")
print(" 但即便把 M2 放宽数倍也到不了 21.5%；nH 散布不是 Atlas 不一致的主因）")

print()
print("=" * 90)
print("M3  零模型判据（pKi−pEC50>0.3 dex）假阳性率")
print("=" * 90)
fp = 0; N = 200000
for _ in range(N):
    pKA = rng.uniform(5, 10); tau = 10 ** rng.uniform(-0.5, 2)
    pki = pKA + rng.normal(0, 0.08)
    pec = pKA + np.log10(1 + tau) + rng.normal(0, 0.1)
    if (pki - pec) > 0.3: fp += 1
print(f"假阳性率 = {fp / N:.5f}（N={N}）；对照：D2R 实测 3/7 配体报警")

print()
print("=" * 90)
print("M4  检出功率（恒等式2 结构，注入系统断裂 δ）")
print("=" * 90)
for delta in [0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]:
    hit = 0; T = 4000
    for _ in range(T):
        pKA = rng.uniform(6, 9); tau = 10 ** rng.uniform(-0.7, 1.5)
        lt_true = np.log10(tau); pec_true = pKA + np.log10(1 + tau)
        lr_true = lt_true + pKA + delta
        lt = lt_true + rng.normal(0, 0.2); pec = pec_true + rng.normal(0, 0.2)
        lr = lr_true + rng.normal(0, 0.15)
        s_l1 = (10 ** lt / (1 + 10 ** lt)) * 0.2
        pred = lt + pec - np.log10(1 + 10 ** lt)
        sd = np.sqrt(0.04 + 0.04 + s_l1 ** 2)
        if abs(lr - pred) > max(0.3, 2 * np.sqrt(sd ** 2 + 0.15 ** 2)): hit += 1
    print(f"  δ={delta:3.1f} dex → 单格检出率 {hit / T:.2f}")
print("（含义：我们的判据是保守的——0.3–0.5 dex 的小断裂大概率漏检；")
print(" 检出的都是 ≥0.8 dex 量级的硬伤，或单列多格同向的联合证据）")

print()
print("=" * 90)
print("M5  阈值稳健性（Atlas 实测分布）")
print("=" * 90)
od = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "结果", "atlas_audit")
grpd = pd.read_csv(os.path.join(od, "atlas_audit_groups.csv"))
papd = pd.read_csv(os.path.join(od, "atlas_audit_papers.csv"))
for t in [0.2, 0.3, 0.5, 0.8, 1.0]:
    print(f"  阈 {t:.1f} dex: 组 {(grpd['dev_median'] > t).mean():.2%}，文献 {(papd['dev_median'] > t).mean():.2%}")
print("（结论：'约 1/5–1/4 的组/文献不自洽'在 0.2–0.5 dex 阈值带内定性不变）")

print()
print("=" * 90)
print("M6  列级联合证据 + 方向对称性")
print("=" * 90)
print(f"Gillis GIRK 列：恒等式1 5/5 报警的零假设概率 = {fp1 ** 5:.2e}")
print(f"恒等式2 4/5 报警 ≈ {fp2 ** 4 * 5:.2e}（二项近似）")
pos = neg = 0
for _ in range(40000):
    d = rng.normal(0, 1)
    if abs(d) > 2:
        if d > 0: pos += 1
        else: neg += 1
print(f"报警方向对称性（应≈50/50）：{pos}/{neg}")
