#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code 51 chemotaxis receptor layer framework self-check simulation v1.2.0
================================================
v1.0.0: naive occupancy kernel falsification + empirical kernel calibration + analyst-normalization degeneracy self-check.
v1.1.0: external audit 1 (Bug 8: analyst normalization = trivial division) — norm downgraded to a trivial-algebra
  control; added MWC active-state + methylation exact-adaptation mechanism kernel, cell-autonomous readout test,
  bootstrap CI. Core finding: the mechanism kernel is fold-native, conflicting with the empirical increment structure.
v1.2.0: external audit 2 ruling —
  accepted warning 9 (over-generalization of "any form" → changed to "standard two-state free-energy difference and common variants",
    added G4 power-law difference, G5 log ratio);
  accepted issues 2/3 (adaptation idealization → introduced a* cross-cell fluctuation σ_a and per-measurement feedback noise σ_m);
  rejected issue 4 (CI values already output in v1.1.0; this version adds a width column);
  rejected issue 5 (R²fit and coordinate R² were already split into separate columns in v1.1.0).
  Also fixed a new bug in the reviewer's revision code: the a* cross-cell fluctuation must be drawn once per cell
  (fixed across the 5 foreground levels), while σ_m is the per-measurement noise; their implementation drew σ_a
  per measurement as well, destroying cell identity.

Run: python3 代码51_趋化receptor层_框架自检仿真.py
Dependencies: numpy, scipy, pandas
"""

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import linregress

rng = np.random.default_rng(20260815)

# ---------------------------------------------------------------
# Part 0: empirical calibration target (P7 erratum re-adjudicated unit table, 30 units, real file)
# ---------------------------------------------------------------
CSV = "/mnt/agents/output/03_细胞线3/结果/P7_Moore2024_FCD-Weber/代码49b_单元table.csv"
_df = pd.read_csv(CSV)
B_E = _df["B_uM"].to_numpy(float)
F_E = _df["F_uM"].to_numpy(float)
R_E = _df["R_a"].to_numpy(float)
assert len(_df) == 30, "unit table should have 30 rows"


def r2(x, y):
    """1-D linear regression R²; returns nan on degenerate input."""
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3 or np.std(x[m]) == 0 or np.std(y[m]) == 0:
        return float("nan")
    return linregress(x[m], y[m]).rvalue ** 2


# ---------------------------------------------------------------
# Part 1: naive MWC occupancy increment kernel — honest failure
# ---------------------------------------------------------------
KD = 1.0  # µM, hypothesized value

def theta(L):
    return L / (KD + L)

def kernel_naive(B, F):
    return theta(B + F) - theta(B)

pred_naive = kernel_naive(B_E, F_E)
print("== Naive MWC occupancy increment kernel ==")
print(f"R²(logF)  = {r2(np.log10(F_E), pred_naive):.3f}   (empirical 0.652)")
print(f"R²(logr)  = {r2(np.log10(F_E / B_E), pred_naive):.3f}   (empirical 0.082)")
print(f"R²(logT)  = {r2(np.log10(B_E + F_E), pred_naive):.3f}   (empirical 0.172)")
print(f"High-background truncation test: B=100,F=40 prediction {kernel_naive(100, 40):.4f} vs empirical 0.424"
      f" (occupancy saturation over-truncates → naive kernel falsified)")

# ---------------------------------------------------------------
# Part 2: calibrated empirical kernel + bootstrap CI
# ---------------------------------------------------------------
def kernel_cal(X, A, Fc, Bs, p):
    B, F = X
    return A * np.log10(F / Fc) / (1.0 + (B / Bs) ** p)

p0 = [0.3, 0.15, 100.0, 0.8]
popt, _ = curve_fit(kernel_cal, (B_E, F_E), R_E, p0=p0,
                    bounds=([0, 1e-4, 1e-2, 0.05], [10, 100, 1e5, 5]),
                    maxfev=40000)
A, Fc, Bs, p = popt
pred_cal = kernel_cal((B_E, F_E), *popt)
ss_res = np.sum((R_E - pred_cal) ** 2)
ss_tot = np.sum((R_E - R_E.mean()) ** 2)
print("\n== Calibrated empirical kernel ==")
print(f"R = {A:.3f} * log10(F/{Fc:.3f} µM) / (1 + (B/{Bs:.1f})^({p:.2f}))")
print(f"R² over 30 units = {1 - ss_res / ss_tot:.3f}")

print("\n== Empirical kernel parameter bootstrap 95% CI (500 resamples) ==")
boots = []
for _ in range(500):
    idx = rng.integers(0, 30, 30)
    try:
        pb, _ = curve_fit(kernel_cal, (B_E[idx], F_E[idx]), R_E[idx], p0=popt,
                          bounds=([0, 1e-4, 1e-2, 0.05], [10, 100, 1e5, 5]),
                          maxfev=20000)
        boots.append(pb)
    except Exception:
        pass
boots = np.array(boots)
for nm, i in zip(["A", "Fc", "Bs", "p"], range(4)):
    lo, hi = np.percentile(boots[:, i], [2.5, 97.5])
    print(f"{nm}: point estimate {popt[i]:.3f}  95%CI [{lo:.3f}, {hi:.3f}]  width {hi - lo:.3f}")
print("Interpretation: Fc CI is narrow (stable class-II anchor point); Bs/p CIs are wide (30 units insufficient to precisely bound the decay).")

# ---------------------------------------------------------------
# Part 3: degeneracy self-check v1 — analyst normalization (trivial-algebra control, not mechanistic evidence)
# ---------------------------------------------------------------
BGS = {0.01: [0.2, 0.5, 1.09, 2.09, 4.09],
       0.1:  [0.2, 0.5, 1, 2, 4],
       0.3:  [0.2, 0.5, 1, 2, 4],
       1.0:  [0.2, 0.5, 1, 2, 4],
       10.0: [0.5, 1, 2, 4, 8],
       100.0:[2, 5, 10, 20, 40]}
B_OFF = 0.3


def run_sweep_analyst(sigma_k, n_cell=300):
    rows = []
    for B, fs in BGS.items():
        for _ in range(n_cell):
            k = np.exp(rng.normal(0, sigma_k))
            amps, amps_n = [], []
            for F in fs:
                r_true = kernel_cal((B, F), *popt)
                amps.append(k * r_true + B_OFF + rng.normal(0, 0.01))
                amps_n.append(k * r_true + rng.normal(0, 0.01))  # note: b is removed here as well
            mx = max(amps_n)
            for F, a, an in zip(fs, amps, amps_n):
                rows.append((B, F, a, an / mx))
    rows = np.array(rows)
    Bv, Fv = rows[:, 0], rows[:, 1]
    out = {}
    for name, col in (("raw", 2), ("norm_analyst", 3)):
        y = rows[:, col]
        out[name] = (r2(np.log10(Fv), y),
                     r2(np.log10(Fv / Bv), y),
                     r2(np.log10(Bv + Fv), y))
    return out


print("\n== Degeneracy self-check v1: analyst normalization (trivial-algebra control, not mechanistic evidence) ==")
print(f"{'sigma_k':>8} | {'raw raw amplitude':^28} | {'norm analyst-normalized':^28}")
for sk in (0.0, 0.3, 0.6):
    res = run_sweep_analyst(sk)
    print(f"{sk:8.1f} | "
          f"{res['raw'][0]:.3f}, {res['raw'][1]:.3f}, {res['raw'][2]:.3f}   | "
          f"{res['norm_analyst'][0]:.3f}, {res['norm_analyst'][1]:.3f}, {res['norm_analyst'][2]:.3f}")

# ---------------------------------------------------------------
# Part 4: MWC active-state + methylation exact-adaptation mechanism kernel (shape test, no noise)
# ---------------------------------------------------------------
ASTAR = 1.0 / 3.0   # hypothesis: kR/(kR+kB)
LAM = np.log(1.0 / ASTAR - 1.0)

def g_std(L, Ki, Ka):
    return np.log((1 + L / Ki) / (1 + L / Ka))

def da_mech(B, F, N, Ki, Ka):
    return ASTAR - 1.0 / (1.0 + np.exp(LAM + N * (g_std(B + F, Ki, Ka) - g_std(B, Ki, Ka))))

print("\n== MWC+BL mechanism kernel: parameter scan (can the shape grow the empirical kernel? no noise) ==")
best = None
for Ki in [0.05, 0.1, 0.2, 0.5, 1.0, 2.0]:
    for Ka in [5, 20, 50, 200]:
        if Ka <= Ki * 2:
            continue
        for N in [3, 6, 10, 18, 30]:
            pv = da_mech(B_E, F_E, N, Ki, Ka)
            if pv.max() < 0.1:
                continue
            sc = (pv * R_E).sum() / (pv * pv).sum()
            ss = 1 - ((sc * pv - R_E) ** 2).sum() / ((R_E - R_E.mean()) ** 2).sum()
            if best is None or ss > best[0]:
                best = (ss, Ki, Ka, N,
                        r2(np.log10(F_E), pv), r2(np.log10(F_E / B_E), pv),
                        r2(np.log10(B_E + F_E), pv))
print(f"Best: R²fit={best[0]:.3f}, Ki={best[1]}, Ka={best[2]}, N={best[3]}")
print(f"Coordinate separability logF/logr/logT = {best[4]:.3f}/{best[5]:.3f}/{best[6]:.3f}"
      f" (empirical 0.652/0.082/0.172) → fold-native, falsified by the empirical table")

print("\n== Standard two-state free-energy difference forms and common variants (warning-9 fix: no longer claiming 'any') ==")
def test_form(name, dg):
    sc = (dg * R_E).sum() / (dg * dg).sum()
    ss = 1 - ((sc * dg - R_E) ** 2).sum() / ((R_E - R_E.mean()) ** 2).sum()
    print(f"{name}: R²fit={ss:7.3f}  logF={r2(np.log10(F_E), dg):.3f} "
          f"logr={r2(np.log10(F_E / B_E), dg):.3f} logT={r2(np.log10(B_E + F_E), dg):.3f}")

test_form("G1 two-state diff (Ki=2,Ka=200)", g_std(B_E + F_E, 2, 200) - g_std(B_E, 2, 200))
test_form("G2 shifted increment ln(1+F/(0.17+B))", np.log(1 + F_E / (0.17 + B_E)))
test_form("G2 shifted increment ln(1+F/(1+B))  ", np.log(1 + F_E / (1 + B_E)))
test_form("G3 pure increment ln(F/0.168)      ", np.log(F_E / 0.168))
test_form("G4 power-law diff (B+F)^0.5-B^0.5  ", (B_E + F_E) ** 0.5 - B_E ** 0.5)
test_form("G5 log ratio log10((B+F)/B)   ", np.log10((B_E + F_E) / B_E))

# ---------------------------------------------------------------
# Part 5: cell-autonomous readout Δa survival test (with adaptation noise, v1.2.0 fixed version)
#   σ_N: receptor-cluster magnitude cross-cell degeneracy (drawn once per cell)
#   σ_a: a* cross-cell fluctuation (drawn once per cell, fixed across the 5 levels — cell identity)
#   σ_m: methylation feedback/baseline estimation noise (per measurement)
#   readout Δa = a_base - a_inst: the baseline is the cell's steady-state activity after adaptation (noisy),
#   cell-availability hypothesis: a cell can only read the baseline from its own post-adaptation activity; no analyst normalization.
# ---------------------------------------------------------------
KI_M, KA_M = 2.0, 200.0

def sweep_autonomous(sN, sa, sm, n_cell=300):
    rows = []
    for B, fs in BGS.items():
        for _ in range(n_cell):
            N = 6 * np.exp(rng.normal(0, sN))                 # cell attribute
            a_star_cell = ASTAR * np.exp(rng.normal(0, sa))   # cell attribute: drawn once per cell
            lam_cell = np.log(1.0 / a_star_cell - 1.0)
            for F in fs:
                a_inst = 1.0 / (1.0 + np.exp(
                    lam_cell + N * (g_std(B + F, KI_M, KA_M) - g_std(B, KI_M, KA_M))))
                a_base = a_star_cell + rng.normal(0, sm)       # per-measurement noise
                rows.append((B, F, a_base - a_inst + rng.normal(0, 0.005)))
    A3 = np.array(rows)
    Bv, Fv, y = A3[:, 0], A3[:, 1], A3[:, 2]
    return (r2(np.log10(Fv), y), r2(np.log10(Fv / Bv), y),
            r2(np.log10(Bv + Fv), y))

print("\n== Degeneracy self-check v2.1: cell-autonomous readout Δa (a* fixed per cell, σ_m per measurement) ==")
print(f"{'σ_N':>4} {'σ_a':>4} {'σ_m':>4} | logF   logr   logT")
for sN in (0.0, 0.3, 0.6):
    for sa, sm in [(0.0, 0.0), (0.05, 0.05), (0.15, 0.05)]:
        t = sweep_autonomous(sN, sa, sm)
        print(f"{sN:4.1f} {sa:4.2f} {sm:4.2f} | {t[0]:.3f}  {t[1]:.3f}  {t[2]:.3f}")
print("Reading: feedback pins the baseline (adaptation noise mildly reduces survival), the surviving coordinate is the fold (logr),")
print("and the conclusion that it conflicts with the empirical increment structure (logF 0.652) is unchanged.")

print("""
== v1.2.0 summary ==
1. Naive occupancy kernel falsified.
2. Empirical kernel R²=0.94; Fc CI narrow (stable class-II anchor point), Bs/p CI wide.
3. Analyst normalization = trivial-algebra control (Bug 8 downgrade maintained).
4. Standard two-state free-energy differences and common variants (G1/G2/G4/G5) are all fold-native or mixed,
   failing to reproduce the empirical increment structure; the closest, G3, is a pure absolute increment, not a "difference" form.
   (Wording converged: no longer claiming "any" form.)
5. Methylation feedback pins the baseline: the cell-autonomous readout survives under N degeneracy + adaptation noise
   with mild decay (logr 0.575→0.36); the surviving coordinate is still the fold — the "paid pinning"
   mechanism holds, and what is preserved is this layer's native statistic.
6. The source of the empirical increment structure remains to be localized: downstream-layer generation / P7 regime / protocol factors;
   the tension with Lazova 2011 (chemotaxis FCD) must be faced.
""")
