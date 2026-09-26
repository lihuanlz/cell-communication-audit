#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code 61 row-wise Hill steepness test: the adjudication experiment for the asymmetric-cooperativity conjecture v1.0.0
================================================
Task: adjudicate the "asymmetric cooperativity" conjecture — if T-state cooperativity is weak (N≈2-4, explaining the
K1/2 Weber line) while R-state cooperativity is strong (N≈12-39, explaining the amplitude table), then the apparent
Hill coefficient of the measured dose curves should exceed the range any single-N MWC model can produce in the same window.

Design (replacing the user's original "bin by amplitude, measure local slope" plan — that plan is contaminated by the
geometric artifact that saturation curves are naturally flat at both ends, see §22):
  §1 Per background level, Hill-fit the whole-population median dose curve da(T)=dmax·Tⁿ/(Kⁿ+Tⁿ):
     (a) dmax free (three parameters); (b) dmax=1 fixed (saturated-pulse calibration meaning).
  §2 Model-side control: standard single-N MWC with literature parameters (Ki=0.81,Ka=200,N=6) and with
     amplitude-optimal parameters (Ki=30.7,Ka=133,N=39), Hill-fitted in exactly the same concentration windows,
     giving the apparent-n envelope of single-N models.
  Criterion: if measured n(B) systematically exceeds the envelope -> asymmetric cooperativity is supported;
        if it falls inside the envelope -> falsified, this branch is sealed.

Results (v1.0.0 live run):
  Measured n(dmax free): 0.91 / 1.40 / 1.67 / 2.34 / 3.14 / 13.63 / 16.85
                    (B = 0 / 0.01 / 0.1 / 0.3 / 1 / 10 / 100)
  Measured n(dmax=1):   0.69 / 0.79 / 0.87 / 1.09 / 1.29 / 3.12 / 2.11
                    (high-B rows R²=0.81/0.69; fixed-asymptote Hill fits poorly at high B)
  Literature model n:   1.96 / 2.01 / 2.47 / 3.35 / 5.78 / 11.76 / 12.10
  Amplitude-optimal n:  1.25 / 1.28 / 1.51 / 1.98 / 3.52 / 14.28 / 26.33
  -> The rise of measured n(B) (~1->17) has the same shape as the rise of the single-N MWC same-window apparent n,
    and falls everywhere inside or below the envelope of the two single-N models.
  -> Asymmetric-cooperativity conjecture falsified: no second N is needed. This branch is sealed.

Note: model dmax_app->0 at high B is because the window covers only the bottom of the incremental curve (the same-source
phenomenon as the §21 tail-statistic warning); it does not affect the envelope comparison of n.

Data: Moore 2024 Dryad doi:10.5061/dryad.nvx0k6dzz (CC0)
Run: python3 code61_rowwise_Hill_steepness_test.py
Depends: numpy, scipy, code54 (same directory)
"""

import os
import importlib.util
import warnings
from collections import defaultdict

import numpy as np
import scipy.io as sio
from scipy.optimize import curve_fit

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "c54", os.path.join(HERE, "code54_K12_joint_fit_rupture_quantification.py"))
c54 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c54)

ROOT = c54.ROOT
FILES = c54.FILES

ASTAR = 1.0 / 3.0
LAM = np.log(1.0 / ASTAR - 1.0)

def gL(L, Ki, Ka):
    L = np.maximum(L, 0.0)
    return np.log((1 + L / Ki) / (1 + L / Ka))

def a_post(L, Ki, Ka, N):
    return 1.0 / (1.0 + np.exp(LAM + N * gL(L, Ki, Ka)))

def mwc_da(Ts, B, Ki, Ka, N):
    """Incremental curve of the single-N standard MWC (amax=1, a*=1/3 fixed)."""
    return a_post(B, Ki, Ka, N) - a_post(np.asarray(Ts), Ki, Ka, N)

def hill(T, K, n, dmax):
    T = np.asarray(T)
    return dmax * T**n / (K**n + T**n)

def hill_fix1(T, K, n):
    return hill(T, K, n, 1.0)

# ---------------------------------------------------------------
# Data: per-cell (T=B+F, da) pooling (F=sv−B, consistent with the code-59 fixed version)
# ---------------------------------------------------------------
def pooled_curves():
    out = defaultdict(list)
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
                F = sv - B                    # same-ligand files: s includes background
                if F <= 0:
                    continue
                lev[round(F, 4)].append(float(np.median(a[row, pre[-4:]]))
                                        - float(np.median(a[row, stim[-6:]])))
            good = {k: float(np.median(v)) for k, v in lev.items() if len(v) >= 4}
            if len(good) >= 3:
                for F, v in good.items():
                    out[B].append((B + F, v))
    return out

def fit_hill(Ts, das, fix_dmax1=False):
    if fix_dmax1:
        popt, _ = curve_fit(hill_fix1, Ts, das,
                            p0=[np.median(Ts), 2.0], maxfev=20000)
        K, n, dmax = popt[0], popt[1], 1.0
        pred = hill_fix1(Ts, K, n)
    else:
        popt, _ = curve_fit(hill, Ts, das,
                            p0=[np.median(Ts), 2.0, max(das)], maxfev=20000)
        K, n, dmax = popt
        pred = hill(Ts, K, n, dmax)
    ss = 1 - np.sum((das - pred)**2) / np.sum((das - das.mean())**2)
    return K, n, dmax, ss

def main():
    pooled = pooled_curves()

    print("§1 measured row-wise Hill fits (whole-population median curves, T=B+F axis)")
    print(f"{'B':>7} {'n(free)':>8} {'K(free)':>9} {'dmax':>6} {'R²':>6}"
          f" | {'n(dmax=1)':>10} {'K(dmax=1)':>10} {'R²':>6}")
    rows = {}
    for B in sorted(pooled):
        pts = np.array(pooled[B])
        Tlev = sorted(set(np.round(pts[:, 0], 4)))
        med = np.array([np.median(pts[np.round(pts[:, 0], 4) == t, 1]) for t in Tlev])
        Tlev = np.array(Tlev)
        rows[B] = (Tlev, med)
        K1, n1, d1, r1 = fit_hill(Tlev, med)
        K2, n2, _, r2 = fit_hill(Tlev, med, fix_dmax1=True)
        print(f"{B:7.2f} {n1:8.2f} {K1:9.2f} {d1:6.3f} {r1:6.3f}"
              f" | {n2:10.2f} {K2:10.2f} {r2:6.3f}")

    print("\n§2 model-side same-window apparent Hill coefficients (dmax free)")
    print(f"{'B':>7} {'lit N=6':>9} {'amp-opt N=39':>13}")
    for B in sorted(rows):
        Ts, _ = rows[B]
        line = f"{B:7.2f}"
        for Ki, Ka, N in [(0.81, 200.0, 6.0), (30.7, 133.0, 39.0)]:
            Dm = mwc_da(Ts, B, Ki, Ka, N)
            popt, _ = curve_fit(hill, Ts, Dm,
                                p0=[np.median(Ts), 2.0, max(Dm)], maxfev=20000)
            line += f" {popt[1]:9.2f}" if N == 6 else f" {popt[1]:13.2f}"
        print(line)

if __name__ == "__main__":
    main()
