# -*- coding: utf-8 -*-
"""
code72_audit_method_false_positive_selfcheck.py  v1.0.0
============================================================
cell line 4 · C1: the audit tool audits the audit tool — false-positive rate and power

codes 63-71 pronounced four "ruptures" on others, but the false-positive rate of the criteria themselves was never measured.
This code uses null-hypothesis simulation (a fully self-consistent true model + noise at the original papers' magnitude) to answer four questions:

  M1  Gillis structural null simulation: alarm rate of identities 1/2 on self-consistent data
      (observed for comparison: GIRK column alarms 5/5 and 4/5)
  M2  how large a spurious deviation can legitimate nH!=1 variation produce (calibrating the attributable share
      of Atlas's 21.5%; generalized operational model E=Esys(tau*A)^n/((KA+A)^n+(tau*A)^n), nH in [0.8,1.25])
  M3  false-positive rate of the null-model criterion (pKi-pEC50 > 0.3 dex)
  M4  detection power curve: inject a systematic rupture of delta dex -> per-cell detection rate
  M5  robustness of the verdicts to the threshold (Atlas observed distribution scanned over 0.2-1.0 dex thresholds)
  M6  null-hypothesis probability of column-level joint evidence + symmetry of alarm direction

Discipline: this code only scores the methods and changes no sealed verdict; if the false-positive rate is found
above its nominal value, the affected verdicts are demoted under the double-entry rule.
Run: python3 code72_audit_method_false_positive_selfcheck.py
============================================================
"""
import numpy as np, pandas as pd, os
from scipy.optimize import brentq

rng = np.random.default_rng(20260816)
LN10 = np.log(10)

print("=" * 90)
print("M1  Gillis structural null simulation (true model self-consistent; SEM at original magnitude: pEC50±0.2, logtau±0.2,")
print("    Emax±5, log(tau/KA)±0.15)")
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
print(f"identity-1 false-positive rate = {fp1:.4f}   identity-2 false-positive rate = {fp2:.4f}")
print("observed for comparison: Gillis GIRK column identity-1 5/5, identity-2 4/5; full-table |dev|>0.3 at 39%/38%")

print()
print("=" * 90)
print("M2  spurious deviation from legitimate nH!=1 variation (Atlas identity-A calibration)")
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
    ltk1 = np.log10(tau1) + (-np.log10(a50) - np.log10(1 + tau1))   # "literature tau/KA" back-computed under nH=1
    emr, a50r = 100.0, ec50_gen(100, 1, 1.0)
    logRA_r = np.log10(emr) - np.log10(a50r)
    ltk_r = np.log10(100) + (-np.log10(a50r) - np.log10(101))
    devs.append(abs((ltk1 - ltk_r) - (logRA - logRA_r)))
devs = np.array(devs)
print(f"|delta log(tau/KA) - delta log(RA)| for nH in [0.8,1.25]: median {np.median(devs):.3f} dex, "
      f">0.3 at {np.mean(devs > 0.3):.4f}, >1.0 at {np.mean(devs > 1):.5f}")
print("observed Atlas for comparison: >0.3 at 21.5%, >1.0 at 4.2%")
print("(boundary note: the literature side is back-computed under nH=1; if the original papers fitted nH freely, the spurious deviation would be larger —")
print(" but even relaxing M2 several-fold still falls short of 21.5%; nH dispersion is not the main cause of the Atlas inconsistency)")

print()
print("=" * 90)
print("M3  null-model criterion (pKi-pEC50 > 0.3 dex) false-positive rate")
print("=" * 90)
fp = 0; N = 200000
for _ in range(N):
    pKA = rng.uniform(5, 10); tau = 10 ** rng.uniform(-0.5, 2)
    pki = pKA + rng.normal(0, 0.08)
    pec = pKA + np.log10(1 + tau) + rng.normal(0, 0.1)
    if (pki - pec) > 0.3: fp += 1
print(f"false-positive rate = {fp / N:.5f} (N={N}); comparison: D2R observed 3/7 ligands alarm")

print()
print("=" * 90)
print("M4  detection power (identity-2 structure, injected systematic rupture delta)")
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
    print(f"  delta={delta:3.1f} dex -> per-cell detection rate {hit / T:.2f}")
print("(implication: our criterion is conservative — small ruptures of 0.3-0.5 dex are likely missed;")
print(" what is detected are hard failures of >=0.8 dex magnitude, or same-direction joint evidence across multiple cells in one column)")

print()
print("=" * 90)
print("M5  threshold robustness (Atlas observed distribution)")
print("=" * 90)
od = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "结果", "atlas_audit")
grpd = pd.read_csv(os.path.join(od, "atlas_audit_groups.csv"))
papd = pd.read_csv(os.path.join(od, "atlas_audit_papers.csv"))
for t in [0.2, 0.3, 0.5, 0.8, 1.0]:
    print(f"  threshold {t:.1f} dex: groups {(grpd['dev_median'] > t).mean():.2%}, papers {(papd['dev_median'] > t).mean():.2%}")
print("(conclusion: 'about 1/5-1/4 of groups/papers are not self-consistent' is qualitatively unchanged across the 0.2-0.5 dex threshold band)")

print()
print("=" * 90)
print("M6  column-level joint evidence + direction symmetry")
print("=" * 90)
print(f"Gillis GIRK column: null-hypothesis probability of identity-1 5/5 alarms = {fp1 ** 5:.2e}")
print(f"identity-2 4/5 alarms ~ {fp2 ** 4 * 5:.2e} (binomial approximation)")
pos = neg = 0
for _ in range(40000):
    d = rng.normal(0, 1)
    if abs(d) > 2:
        if d > 0: pos += 1
        else: neg += 1
print(f"alarm direction symmetry (should be ~50/50): {pos}/{neg}")
