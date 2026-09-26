#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code 56 FRET readout-layer nonlinearity test v1.0.0
================================================
Candidate 2: does the coordinate transform happen at the readout layer? Upstream fixed to standard MWC + precise adaptation
(Ki, Ka, N free), CheY-P linear in activity; three CheY-P->FRET mappings tested:
  V1 saturating-binding map S=Y/(1+Y), Y=c·a (hyperbola of CheY-P/CheZ binding), absolute response;
  V2 same map + pulse normalized by per-background saturation (matching the calibration of the a field);
  V3 power-law map S∝Y^q (q>1 expansive / q<1 compressive), testing the direction of the coordinate transform.
Discrimination as before: amplitude R²>0.9 and K1/2 error<0.1 dex. K1/2 is the absolute half-amplitude da=0.5;
model-prediction censoring penalized 1 dex. least_squares best of 14 starts.

Live-run results (v1.0.0):
  V1 saturating map: c->0.34-0.37 (CheY-P always sits in the near-linear readout region; saturation is not
    demanded by the data), amax hits boundary 5. lam=1: R²=0.946/err=0.256.
    -> the saturating direction (compression) is opposite to what is needed: the data are "more incremental"
      than the standard model; compression only makes logF more like logF.
  V2 background-normalized: worse (lam=1: R²=0.754/err=0.206, c->0.05 hits lower bound).
  V3 power-law map: q->2.0-2.1 (expansive) appears robustly, lam=1: R²=0.950/err=0.257,
    amax hits boundary 5. Better than standard MWC (0.831/0.224 at the same lam, higher amplitude),
    but not better than M2c imperfect adaptation + single threshold (0.953/0.182), and does not enter the acceptance region.
Verdict: none of the three readout-layer nonlinearity forms can close the fracture. The direction the data need is "expansive"
  (q≈2 robustly selected), but a readout-layer transform alone is insufficient; M2c remains the current best
  (R²=0.953/err=0.182, Ka->inf single-threshold logarithmic sensor).
Significance: candidate 2 is essentially ruled out as the sole main cause. Closing the fracture requires the receptor layer itself
  to have a "low-B plateau K1/2≈2.5 µM + steep rise at B=10 to 14 µM" structure — pointing to candidate 1
  (effective parameters changing with methylation state/background) or a combination of readout expansion and receptor logarithmization.

Data: Moore 2024 Dryad doi:10.5061/dryad.nvx0k6dzz (CC0)
Run: python3 code56_FRET_readout_nonlinearity_test.py
Depends: numpy, scipy, pandas; code54 in the same directory (data extraction reused)
"""

import importlib.util
import os
import warnings

import numpy as np
from scipy.optimize import brentq, least_squares

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "c54", os.path.join(HERE, "code54_K12_joint_fit_rupture_quantification.py"))
c54 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c54)

ASTAR = 1.0 / 3.0
LAM = np.log(1.0 / ASTAR - 1.0)

def gL(L, Ki, Ka):
    return np.log((1 + L / Ki) / (1 + L / Ka))

def load_data():
    k12, cens, tot = c54.extract_K12()
    K12_B = np.array(sorted(set(list(k12) + list(cens))))
    K12_obs = np.array([10 ** np.nanmedian(np.log10(np.array(k12.get(B, [np.nan]))))
                        for B in K12_B])
    import pandas as pd
    df = pd.read_csv(c54.CSV)
    return (df["B_uM"].to_numpy(float), df["F_uM"].to_numpy(float),
            df["R_a"].to_numpy(float), K12_B, K12_obs)

def a_post_of(B, F, Ki, Ka, N):
    return 1 / (1 + np.exp(LAM + N * (gL(B + F, Ki, Ka) - gL(B, Ki, Ka))))

def da_sat(B, F, Ki, Ka, N, c, amax):          # V1
    y0, y1 = c * ASTAR, c * a_post_of(B, F, Ki, Ka, N)
    return amax * (y0 / (1 + y0) - y1 / (1 + y1))

def da_satn(B, F, Ki, Ka, N, c, amax):         # V2 (amax placeholder, unused)
    y0, y1 = c * ASTAR, c * a_post_of(B, F, Ki, Ka, N)
    sat = y0 / (1 + y0)
    return (y0 / (1 + y0) - y1 / (1 + y1)) / sat if sat > 1e-12 else 0.0

def da_pow(B, F, Ki, Ka, N, q, amax):          # V3
    return amax * (ASTAR ** q - a_post_of(B, F, Ki, Ka, N) ** q)

def fit(da_fn, B_E, F_E, R_E, K12_B, K12_obs, lam, nstart=14):
    def K12(B, p):
        if da_fn(B, 1e6, *p) < 0.5:
            return np.nan
        try:
            return brentq(lambda F: da_fn(B, F, *p) - 0.5, 1e-6, 1e6, xtol=1e-8)
        except Exception:
            return np.nan
    def resid(t):
        p = tuple(np.exp(t))
        ra = (da_fn(B_E, F_E, *p) - R_E) / np.std(R_E)
        rk = [np.log10(K12(B, p) / o) if np.isfinite(K12(B, p)) else 1.0
              for B, o in zip(K12_B, K12_obs)]
        return np.concatenate([ra, np.sqrt(lam) * np.array(rk) / 0.5])
    best = None
    lo = np.log([0.05, 5, 1, 0.05, 0.3]); hi = np.log([300, 1e6, 60, 300, 5])
    for seed in range(nstart):
        r_ = np.random.default_rng(seed)
        xs = np.log([2, 200, 6, 3, 1.0]) + r_.normal(0, 0.8, 5)
        try:
            sol = least_squares(resid, xs, bounds=(lo, hi), max_nfev=15000)
            ss = float(np.sum(resid(sol.x) ** 2))
            if best is None or ss < best[0]:
                best = (ss, sol.x)
        except Exception:
            pass
    p = tuple(np.exp(best[1]))
    pr = da_fn(B_E, F_E, *p)
    r2a = 1 - ((pr - R_E) ** 2).sum() / ((R_E - R_E.mean()) ** 2).sum()
    errs = [abs(np.log10(K12(B, p) / o)) if np.isfinite(K12(B, p)) else 1.0
            for B, o in zip(K12_B, K12_obs)]
    return r2a, float(np.mean(errs)), p

def main():
    B_E, F_E, R_E, K12_B, K12_obs = load_data()
    for name, fn, qmode in [("V1 saturating map", da_sat, False),
                            ("V2 background-normalized", da_satn, False),
                            ("V3 power-law map", da_pow, True)]:
        print(f"== {name} ==")
        lo_q = 0.2 if qmode else 0.05
        for lam in [0.0, 1.0, 10.0]:
            if qmode:
                # V3 parameter order (Ki,Ka,N,q,amax)
                def resid_q(t):
                    p = tuple(np.exp(t))
                    ra = (da_pow(B_E, F_E, *p) - R_E) / np.std(R_E)
                    def K12(B, p):
                        if da_pow(B, 1e6, *p) < 0.5: return np.nan
                        try: return brentq(lambda F: da_pow(B, F, *p) - 0.5, 1e-6, 1e6, xtol=1e-8)
                        except Exception: return np.nan
                    rk = [np.log10(K12(B, p) / o) if np.isfinite(K12(B, p)) else 1.0
                          for B, o in zip(K12_B, K12_obs)]
                    return np.concatenate([ra, np.sqrt(lam) * np.array(rk) / 0.5])
                best = None
                for seed in range(14):
                    r_ = np.random.default_rng(seed)
                    xs = np.log([2, 200, 6, 1.0, 1.0]) + r_.normal(0, 0.8, 5)
                    try:
                        sol = least_squares(resid_q, xs,
                            bounds=(np.log([0.05, 5, 1, 0.2, 0.3]), np.log([300, 1e6, 60, 4, 5])),
                            max_nfev=15000)
                        ss = float(np.sum(resid_q(sol.x) ** 2))
                        if best is None or ss < best[0]: best = (ss, sol.x)
                    except Exception: pass
                p = tuple(np.exp(best[1]))
                pr = da_pow(B_E, F_E, *p)
                r2a = 1 - ((pr - R_E) ** 2).sum() / ((R_E - R_E.mean()) ** 2).sum()
                def K12p(B):
                    if da_pow(B, 1e6, *p) < 0.5: return np.nan
                    try: return brentq(lambda F: da_pow(B, F, *p) - 0.5, 1e-6, 1e6, xtol=1e-8)
                    except Exception: return np.nan
                errs = [abs(np.log10(K12p(B) / o)) if np.isfinite(K12p(B)) else 1.0
                        for B, o in zip(K12_B, K12_obs)]
                print(f"  λ={lam:5.1f}: R²={r2a:.3f}, K1/2err={np.mean(errs):.3f} | "
                      f"Ki={p[0]:.2f}, Ka={p[1]:.0f}, N={p[2]:.1f}, q={p[3]:.2f}, amax={p[4]:.2f}")
            else:
                r2a, er, p = fit(fn, B_E, F_E, R_E, K12_B, K12_obs, lam)
                print(f"  λ={lam:5.1f}: R²={r2a:.3f}, K1/2err={er:.3f} | "
                      f"Ki={p[0]:.2f}, Ka={p[1]:.0f}, N={p[2]:.1f}, c={p[3]:.3f}, amax={p[4]:.2f}")

    print("""
== Summary ==
V1/V2 saturating maps: rejected by the data (c->small, near-linear region; compression direction opposite to need).
V3 power-law map: q≈2 (expansive) appears robustly, reaching R²=0.950/err=0.257 at lam=1,
  better than plain standard MWC but not better than M2c (0.953/0.182); does not enter the acceptance region.
Readout-layer nonlinearity as the sole main cause is essentially ruled out; the data need a transform in the "expansive" direction;
fracture closure points to receptor-layer structure (candidate 1 gain retuning) or a combination with readout expansion.
""")

if __name__ == "__main__":
    main()
