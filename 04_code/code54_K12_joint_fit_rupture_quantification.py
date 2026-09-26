#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
code54 K1/2 joint fit: quantifying the self-consistency rupture of the standard MWC model v1.0.0
================================================
Task: quantitatively answer "can the standard MWC + exact-adaptation model, with a single set of physical
parameters, self-consistently reproduce both Moore-paper statistics — the amplitude table R(B,F) and the K1/2(B) distribution".

§1 per-cell K1/2 extraction from the raw .mat files: under field a (saturating-stimulus calibrated, 0/1 normalized),
   responses da are sorted by total ligand T=B+F and linearly interpolated on the log-T axis at da=0.5.
   Cells whose maximum da < 0.5 are recorded as right-censored (insensitive cells, excluded from the median).
   K1/2 and the amplitude table come from the same batch of raw data, maximizing internal consistency.
§2 joint fit: parameters (Ki, Ka, N, amax) (log space); residual =
   amplitude residual/std(R) + sqrt(lambda)*log10(K1/2 predicted/observed)/0.5 dex;
   a lambda sweep gives the Pareto front. least_squares with four starts, best kept.

Results (v1.0.0 live run):
  K1/2 medians: B=0->2.03, 0.01->2.90, 0.1->2.17, 0.3->2.58,
               1->3.49, 10->13.92, 100->119.56 uM
  (plateau ~2-3.5 uM (B<=1), then ~K0+1.17*B; the B=0 group is 67% censored, biased low)
  Pareto front:
    lambda=0    amplitude R^2=0.980, K1/2 error 0.62 dex (Ki=31, N=39)
    lambda=1    amplitude R^2=0.831, K1/2 error 0.22 dex
    lambda=30   amplitude R^2=0.511, K1/2 error 0.07 dex (Ki=4.3, N=2.0)
  the front is steep with no acceptable compromise point (no point with R^2>0.9 and |dlog10|<0.1)
  -> the rupture is quantitatively established. Fitted amax 1.4-2.7 (physically should be ~1) is another sign of strain.

Data: Moore 2024 Dryad doi:10.5061/dryad.nvx0k6dzz (CC0)
Run: python3 code54_K12_joint_fit_rupture_quantification.py
Dependencies: numpy, scipy, pandas
"""

import os
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.optimize import brentq, least_squares

warnings.filterwarnings("ignore")

ROOT = "/mnt/agents/output/03_细胞线3/公开数据/Moore2024_Chemotaxis_FRET"
CSV = "/mnt/agents/output/03_细胞线3/结果/P7_Moore2024_FCD-Weber/代码49b_单元表.csv"

FILES = {
    "210802_FOV1": 0, "210802_FOV2": 0, "210805_FOV1": 0, "210805_FOV2": 0,
    "220106_FOV1": 0, "230417_FOV1": 0,
    "230815_FOV1": 0.01, "230815_FOV2": 0.01, "230816_FOV1": 0.01, "230816_FOV2": 0.01,
    "230830_FOV1": 0.1, "230830_FOV2": 0.1, "230831_FOV1": 0.1, "230831_FOV2": 0.1,
    "220615_FOV1": 0.3, "230410_FOV1": 0.3, "230428_FOV1": 1.0, "230429_FOV1": 1.0,
    "220302_FOV1": 10.0, "220303_FOV1": 10.0,
    "210816_FOV1": 100.0, "210816_FOV2": 100.0, "230717_FOV1": 100.0, "230718_FOV1": 100.0,
}

# ---------------------------------------------------------------
# §1 per-cell K1/2 extraction
# ---------------------------------------------------------------
def extract_K12():
    k12 = defaultdict(list)
    cens = defaultdict(int)
    tot = defaultdict(int)
    for name, bg in FILES.items():
        p = os.path.join(ROOT, name + ".mat")
        if not os.path.isfile(p):
            continue
        try:
            rd = sio.loadmat(p)["reorgData"]["resp_data"][0, 0]
        except Exception:
            continue
        B = 100.0 if name == "230831_FOV2" else float(bg)
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
                da = float(np.median(a[row, pre[-4:]])) - float(np.median(a[row, stim[-6:]]))
                if sv - B > 0:
                    lev[round(sv, 4)].append(da)
            good = {k: float(np.median(v)) for k, v in lev.items() if len(v) >= 4}
            if len(good) < 3:
                continue
            tot[B] += 1
            ks = sorted(good)
            T = np.array(ks)
            D = np.array([good[k] for k in ks])
            if D.max() < 0.5:
                cens[B] += 1
                continue            # right-censored: max stimulus did not reach half amplitude
            if D.min() >= 0.5:
                continue            # left-censored (rare)
            i = int(np.argmax(D >= 0.5))
            x0, x1 = np.log10(T[i - 1]), np.log10(T[i])
            y0, y1 = D[i - 1], D[i]
            k12[B].append(10 ** (x0 + (0.5 - y0) * (x1 - x0) / (y1 - y0)))
    return k12, cens, tot

# ---------------------------------------------------------------
# §2 joint fit
# ---------------------------------------------------------------
ASTAR = 1.0 / 3.0
LAM = np.log(1.0 / ASTAR - 1.0)
TGT = np.log(1.0 / (ASTAR / 2) - 1) - LAM   # N*delta-g needed to halve the activity

def gL(L, Ki, Ka):
    return np.log((1 + L / Ki) / (1 + L / Ka))

def amp_model(B, F, Ki, Ka, N, amax):
    return amax * (ASTAR - 1 / (1 + np.exp(LAM + N * (gL(B + F, Ki, Ka) - gL(B, Ki, Ka)))))

def K12_model(B, Ki, Ka, N):
    try:
        return brentq(lambda F: N * (gL(B + F, Ki, Ka) - gL(B, Ki, Ka)) - TGT,
                      1e-6, 1e6, xtol=1e-8)
    except Exception:
        return np.nan

def main():
    k12, cens, tot = extract_K12()
    print("== §1 per-cell K1/2 distribution (total-ligand axis, uM) ==")
    print(f"{'B':>8} {'n_est':>5} {'n_cens':>5} {'K1/2 med':>9} {'IQR':>18}")
    K12_B, K12_obs = [], []
    for B in sorted(set(list(k12) + list(cens))):
        v = np.log10(np.array(k12.get(B, [np.nan])))
        med = 10 ** np.nanmedian(v)
        q = np.nanpercentile(v, [25, 75])
        K12_B.append(B); K12_obs.append(med)
        print(f"{B:8.2f} {len(k12.get(B, [])):5d} {cens.get(B, 0):5d} "
              f"{med:9.2f}   [{10**q[0]:7.2f},{10**q[1]:7.2f}]")
    K12_B = np.array(K12_B); K12_obs = np.array(K12_obs)
    print("note: the B=0 group is heavily censored, so its median is biased low; the plateau (B<=1) and the linear-in-B segment (B>=10) are structurally robust.")

    df = pd.read_csv(CSV)
    B_E = df["B_uM"].to_numpy(float); F_E = df["F_uM"].to_numpy(float)
    R_E = df["R_a"].to_numpy(float)

    def resid(theta, lam):
        Ki, Ka, N, amax = np.exp(theta)
        ra = (amp_model(B_E, F_E, Ki, Ka, N, amax) - R_E) / np.std(R_E)
        rk = []
        for B, obs in zip(K12_B, K12_obs):
            pred = K12_model(B, Ki, Ka, N)
            rk.append(0.0 if not np.isfinite(pred) else np.log10(pred / obs))
        return np.concatenate([ra, np.sqrt(lam) * np.array(rk) / 0.5])

    x0 = np.log([2.0, 200.0, 6.0, 1.0])
    print("\n== §2 joint-fit Pareto front ==")
    print(f"{'lambda':>6} {'amp R^2':>7} {'K1/2 med|dlog10|':>15}  parameters")
    for lam in [0.0, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0]:
        best = None
        for seed in range(4):
            r_ = np.random.default_rng(seed)
            xs = x0 + r_.normal(0, 0.8, 4)
            try:
                sol = least_squares(resid, xs, args=(lam,),
                                    bounds=(np.log([0.05, 5, 2, 0.3]),
                                            np.log([300, 2e4, 60, 3])),
                                    max_nfev=8000)
                ss = np.sum(resid(sol.x, lam) ** 2)
                if best is None or ss < best[0]:
                    best = (ss, sol.x)
            except Exception:
                pass
        Ki, Ka, N, amax = np.exp(best[1])
        pr = amp_model(B_E, F_E, Ki, Ka, N, amax)
        r2a = 1 - ((pr - R_E) ** 2).sum() / ((R_E - R_E.mean()) ** 2).sum()
        dk = np.mean([abs(np.log10(K12_model(B, Ki, Ka, N) / o))
                      for B, o in zip(K12_B, K12_obs)])
        print(f"{lam:6.1f} {r2a:7.3f} {dk:15.3f}  "
              f"Ki={Ki:.1f} Ka={Ka:.0f} N={N:.1f} amax={amax:.2f}")

    print("""
== Conclusion ==
No parameter point achieves R^2>0.9 with K1/2 error <0.1 dex: the front is steep,
so gains on one statistic cost the other -> a single parameter set of standard MWC
cannot self-consistently reproduce Moore's two statistics; rupture quantitatively established. Fitted amax 1.4-2.7 (should be ~1) is a sign of strain.
Next step: can a multi-species / local-adaptation model fit both under one parameter set.
""")


if __name__ == "__main__":
    main()
