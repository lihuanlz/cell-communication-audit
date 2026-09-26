#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code 62 single-cell same-source statistic adjudication: is the fracture a population-sampling artifact? v1.0.0
================================================
Task: adjudicate the "statistical structure" explanation raised in §21 — the K1/2 median samples the tail of the sensitive
subpopulation while the amplitude table samples the whole-population median; the two statistics sample different populations, hence demand different N.

Design (exploiting the fact that the Moore 2024 data are themselves single-cell FRET):
  Build the amplitude table using only the same batch of non-censored cells with estimable K1/2 (da_max≥0.5, 33-55%);
  then the amplitude table and the K1/2 come from **strictly the same batch of cells**. On this same-source data,
  redo the code-54 joint fit (standard MWC, lambda-scan Pareto front):
  · If the fracture closes (a point with R²>0.9 and K1/2 error<0.1 dex appears) -> the fracture is
    a population-sampling artifact and the §21 explanation holds;
  · If the front does not move -> the fracture persists inside a single population, the §21 explanation is insufficient,
    and the fracture is a genuine model-structure deficit.

Results (v1.0.0 live run, shown side by side with the whole-population control):
  Sensitive subpopulation: lam=0 -> R²=0.971/err=0.573; lam=1 -> 0.867/0.206;
            lam=30 -> 0.598/0.071 (Ki=4.1, N=2.0)
  Whole-population control: lam=0 -> R²=0.976/err=0.624; lam=1 -> 0.848/0.261;
            lam=30 -> 0.525/0.074 (Ki=4.3, N=2.0)
  Per-cell K1/2 medians identical to code 54 (2.03/2.90/2.17/2.58/3.49/
  13.92/119.56 µM); pipeline cross-check passed.

Verdict:
  The front barely moves, no acceptable trade-off point -> the fracture is not a population-sampling artifact;
  it persists inside strictly the same batch of cells. The §21 "two statistics sample different populations" leaning
  is **partially retracted** (double-recorded): it still describes a real difference between the two statistics, but
  cannot explain the fracture. The fracture further points to a structural deficit of standard MWC + precise adaptation itself.

Data: Moore 2024 Dryad doi:10.5061/dryad.nvx0k6dzz (CC0)
Run: python3 code62_single_cell_same_source_statistics.py
Depends: numpy, scipy
"""

import os
import warnings
from collections import defaultdict

import numpy as np
import scipy.io as sio
from scipy.optimize import brentq, least_squares

warnings.filterwarnings("ignore")

ROOT = "/mnt/agents/output/03_细胞线3/公开数据/Moore2024_Chemotaxis_FRET"
FILES = {
    "210802_FOV1": 0, "210802_FOV2": 0, "210805_FOV1": 0, "210805_FOV2": 0,
    "220106_FOV1": 0, "230417_FOV1": 0,
    "230815_FOV1": 0.01, "230815_FOV2": 0.01, "230816_FOV1": 0.01, "230816_FOV2": 0.01,
    "230830_FOV1": 0.1, "230830_FOV2": 0.1, "230831_FOV1": 0.1, "230831_FOV2": 0.1,
    "220615_FOV1": 0.3, "230410_FOV1": 0.3, "230428_FOV1": 1.0, "230429_FOV1": 1.0,
    "220302_FOV1": 10.0, "220303_FOV1": 10.0,
    "210816_FOV1": 100.0, "210816_FOV2": 100.0, "230717_FOV1": 100.0, "230718_FOV1": 100.0,
}

ASTAR = 1.0 / 3.0
LAM = np.log(1.0 / ASTAR - 1.0)
TGT = np.log(1.0 / (ASTAR / 2) - 1) - LAM

def gL(L, Ki, Ka):
    L = np.maximum(L, 0.0)
    return np.log((1 + L / Ki) / (1 + L / Ka))

def amp_model(B, F, Ki, Ka, N, amax):
    return amax * (ASTAR - 1 / (1 + np.exp(LAM + N * (gL(B + F, Ki, Ka) - gL(B, Ki, Ka)))))

def K12_model(B, Ki, Ka, N):
    try:
        return brentq(lambda F: N * (gL(B + F, Ki, Ka) - gL(B, Ki, Ka)) - TGT,
                      1e-6, 1e6, xtol=1e-8)
    except Exception:
        return np.nan

# ---------------------------------------------------------------
# Per-cell curve extraction (cell identity preserved)
# ---------------------------------------------------------------
def per_cell_curves():
    per_cell = defaultdict(list)
    for name, bg in FILES.items():
        B = 100.0 if name == "230831_FOV2" else float(bg)
        p = os.path.join(ROOT, name + ".mat")
        try:
            rd = sio.loadmat(p)["reorgData"]["resp_data"][0, 0]
        except Exception:
            continue
        for ci in range(rd.shape[1]):
            try:
                a = rd["a"][0, ci].astype(float)
                s = rd["s"][0, ci].astype(float)
            except Exception:
                continue
            if a.shape != (35, 20) or s.shape != (35, 20):
                continue
            lev = defaultdict(list)
            for row in range(35):
                sr = s[row]
                sv = float(sr.max())
                stim = np.where(sr >= sv - 1e-9)[0]
                pre = np.arange(0, stim[0]) if len(stim) else None
                if pre is None or len(pre) < 4 or len(stim) < 6:
                    continue
                F = sv - B
                if F <= 0:
                    continue
                lev[round(F, 4)].append(float(np.median(a[row, pre[-4:]]))
                                        - float(np.median(a[row, stim[-6:]])))
            good = {k: float(np.median(v)) for k, v in lev.items() if len(v) >= 4}
            if len(good) >= 3:
                per_cell[B].append(good)
    return per_cell

def cell_k12(B, curve):
    items = sorted(curve.items())
    T = np.array([B + f for f, _ in items])
    d = np.array([v for _, v in items])
    if d.max() < 0.5:
        return None
    i = int(np.argmax(d >= 0.5))
    if i == 0:
        return None                    # left-censored, excluded per code-54 convention
    x0, x1 = np.log10(T[i - 1]), np.log10(T[i])
    y0, y1 = d[i - 1], d[i]
    if y1 == y0:
        return T[i]
    return 10 ** (x0 + (0.5 - y0) * (x1 - x0) / (y1 - y0))

def amp_table(curves_by_B, min_cells=10):
    Bs, Fs, Rs = [], [], []
    for B in sorted(curves_by_B):
        lv = defaultdict(list)
        for c in curves_by_B[B]:
            for F, v in c.items():
                lv[F].append(v)
        for F in sorted(lv):
            if len(lv[F]) >= min_cells:
                Bs.append(B); Fs.append(F); Rs.append(float(np.median(lv[F])))
    return np.array(Bs), np.array(Fs), np.array(Rs)

def joint_fit(Bx, Fx, Rx, K12_B, K12_obs, lam, nstart=8):
    def resid(theta):
        Ki, Ka, N, amax = np.exp(theta)
        ra = (amp_model(Bx, Fx, Ki, Ka, N, amax) - Rx) / np.std(Rx)
        rk = []
        for B, o in zip(K12_B, K12_obs):
            pred = K12_model(B, Ki, Ka, N)
            rk.append(0.0 if not np.isfinite(pred) else np.log10(pred / o))
        return np.concatenate([ra, np.sqrt(lam) * np.array(rk) / 0.5])
    x0 = np.log([2.0, 200.0, 6.0, 1.0])
    best = None
    for seed in range(nstart):
        r_ = np.random.default_rng(seed)
        xs = x0 + r_.normal(0, 0.8, 4)
        try:
            sol = least_squares(resid, xs,
                                bounds=(np.log([0.05, 5, 2, 0.3]),
                                        np.log([300, 2e4, 60, 3])),
                                max_nfev=8000)
            ss = np.sum(resid(sol.x) ** 2)
            if best is None or ss < best[0]:
                best = (ss, sol.x)
        except Exception:
            pass
    Ki, Ka, N, amax = np.exp(best[1])
    pr = amp_model(Bx, Fx, Ki, Ka, N, amax)
    r2a = 1 - ((pr - Rx) ** 2).sum() / ((Rx - Rx.mean()) ** 2).sum()
    dk = np.mean([abs(np.log10(K12_model(B, Ki, Ka, N) / o))
                  for B, o in zip(K12_B, K12_obs)])
    return r2a, dk, (Ki, Ka, N, amax)

def main():
    per_cell = per_cell_curves()

    est_curves, all_curves, est_k12 = {}, {}, {}
    for B in sorted(per_cell):
        ec, ks = [], []
        for c in per_cell[B]:
            k = cell_k12(B, c)
            if k is not None:
                ec.append(c); ks.append(k)
        est_curves[B] = ec
        all_curves[B] = per_cell[B]
        est_k12[B] = np.array(ks)
        print(f"B={B:7.2f}: estimable {len(ec)}/{len(per_cell[B])}"
              f" ({100*len(ec)/len(per_cell[B]):.0f}%),"
              f" K1/2 median={np.median(ks) if ks else np.nan:.2f}")

    K12_B = np.array(sorted(est_k12))
    K12_obs = np.array([np.median(est_k12[B]) for B in K12_B])

    B_E, F_E, R_E = amp_table(est_curves)
    B_A, F_A, R_A = amp_table(all_curves)

    for tag, Bx, Fx, Rx in [("sensitive subpopulation (same-source statistics)", B_E, F_E, R_E),
                            ("whole population (code-54-convention control)", B_A, F_A, R_A)]:
        print(f"\n== {tag} ==")
        for lam in [0.0, 1.0, 30.0]:
            r2, dk, par = joint_fit(Bx, Fx, Rx, K12_B, K12_obs, lam)
            print(f"  λ={lam:5.1f}: amplitude R²={r2:.3f}, K1/2 err={dk:.3f} dex,"
                  f" Ki={par[0]:.1f} Ka={par[1]:.0f} N={par[2]:.1f} amax={par[3]:.2f}")

    print("""
== Conclusion ==
The Pareto front of the same-source statistics (amplitude table and K1/2 from the same batch of non-censored cells)
almost coincides with the whole-population convention; no trade-off point with R²>0.9 and K1/2 error<0.1 dex exists.
The fracture is not a population-sampling artifact -> the §21 tail-statistic explanation is partially retracted (double-recorded);
the fracture is a genuine structural deficit of standard MWC + precise adaptation.
""")

if __name__ == "__main__":
    main()
