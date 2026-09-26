#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code 57 TCS depletion robustness check + gain retuning (candidate ①) joint test v1.0.0
================================================
Two parts:
§1 TCS depletion-correction robustness check. All previous receptor-layer models used the zero-depletion approximation
   (g(L) takes total ligand directly). TCS unified equation ξ=p/(1-p)+p/κ, κ=Kd/(nR_T).
   Estimate: densest case 1e9 cells/mL × 1e4 receptors/cell → R_T≈17 nM,
   κ≈Ki/0.017 ≈ 48 ≫ 1 (zero-depletion regime). This code makes no explicit approximation; it solves
   the implicit free-ligand equation L_T = L_free + R_T·L_free/(Kd+L_free) (Kd≈Ki),
   substitutes free ligand into g(L), and rescans the Pareto front.
   Actual run: at λ=1, (R², err) = (0.948, 0.257), identical to the zero-depletion control,
   parameters match digit for digit → depletion formally excluded as a fracture source (consistent with the κ≫1 estimate).

§2 candidate ① gain retuning (effective parameters vary with background/methylation state), three minimal parameterizations:
   G-Ki：Ki_eff(B)=Ki0·(1+B/Bg)^h1
   G-N： N_eff(B)=N0·(1+B/Bg)^h2
   G-both: both together.
   Discrimination unchanged: R²>0.9 and K1/2 error <0.1 dex.

Actual-run results (v1.0.0, block-verified: §1 re-run independently with this file's functions, digit-identical;
§2 variants re-run by the kernel over multiple rounds of 8–16 starts; G-Ki identical across two rounds; G-N/G-both
take the best of all runs — multi-start optimization is sensitive to start count; too few starts fall into an
amplitude-only local optimum (R²=0.980/err=1.00), masking the compromise basin; dual-recorded):
  G-Ki：λ=1 → 0.949/0.254；λ=30 → 0.663/0.064
        (at low λ, h→−2.3 lower bound, Bg→large: degenerates to constant Ki; λ=3 per-B shows
         low plateau ~3.7 vs measured 2.0–2.9, B=10/100 each short by 0.22–0.26 dex)
  G-N ：λ=1 → 0.948/0.247；λ=30 → 0.678/0.052
        (N0→60, h→2.3 both hit bounds: the fit wants gain rising without limit with background,
         physically unreasonable; Ka→9–30 µM far from literature anchors)
  G-both: coincides with G-N (h1→0, Ki falls back to constant) → wasted degrees of freedom.

All-model Pareto comparison (λ=1 compromise point / K1/2 endpoint):
  standard MWC          0.831/0.224 ↔ 0.511/0.074
  M2c incomplete adaptation  0.953/0.182 ↔ 0.843/0.089   ← current best
  FRET power-law readout  0.950/0.257 ↔ 0.758/0.100
  gain retuning G-N     0.948/0.247 ↔ 0.678/0.052

Verdict:
1. depletion: formally excluded (front unmoved digit for digit).
2. minimal gain retuning: cannot close the fracture. It can push the K1/2 end to 0.052 dex
   (best in field), at the cost of amplitude R²=0.68 and parameters hitting non-physical bounds
   (gain rising without limit with background). No model enters the acceptance region.
3. the fracture is robust to all of these perturbations: TCS depletion, two subpopulations, incomplete adaptation,
   readout nonlinearity (compression/expansion), smooth gain retuning.
   Remaining paths: P8 independent-data extrapolation; or directly accept that "a statistical structure
   the standard mechanism family cannot reconcile exists inside the Moore data" as this branch's final verdict, dual-recorded.

Data: Moore 2024 Dryad doi:10.5061/dryad.nvx0k6dzz (CC0)
run: python3 代码57_TCS耗竭稳健性与增益重调检验.py
dependencies: numpy, scipy, pandas; 代码54 in same directory (data extraction reused)
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
RT = 0.017  # µM, densest case

def gL(L, Ki, Ka):
    return np.log((1 + L / Ki) / (1 + L / Ka))

def Lfree(LT, Kd, RT=RT):
    b = Kd + RT - LT
    return (-b + np.sqrt(b * b + 4 * Kd * LT)) / 2

def load_data():
    k12, cens, tot = c54.extract_K12()
    K12_B = np.array(sorted(set(list(k12) + list(cens))))
    K12_obs = np.array([10 ** np.nanmedian(np.log10(np.array(k12.get(B, [np.nan]))))
                        for B in K12_B])
    import pandas as pd
    df = pd.read_csv(c54.CSV)
    return (df["B_uM"].to_numpy(float), df["F_uM"].to_numpy(float),
            df["R_a"].to_numpy(float), K12_B, K12_obs)

# ---------- §1 depletion ----------
def da_dep(B, F, Ki, Ka, N, amax, RT=RT):
    Bf, TF = Lfree(B, Ki, RT), Lfree(B + F, Ki, RT)
    return amax * (ASTAR - 1 / (1 + np.exp(LAM + N * (gL(TF, Ki, Ka) - gL(Bf, Ki, Ka)))))

# ---------- §2 gain retuning ----------
def da_gain(B, F, Ki0, Ka, N0, Bg, h1, h2, amax):
    Kie = Ki0 * (1 + B / Bg) ** h1
    Ne = N0 * (1 + B / Bg) ** h2
    return amax * (ASTAR - 1 / (1 + np.exp(LAM + Ne * (gL(B + F, Kie, Ka) - gL(B, Kie, Ka)))))

def fit_joint(da_fn, pnames, x0, lo, hi, data, lam, nstart=12):
    B_E, F_E, R_E, K12_B, K12_obs = data
    def K12(B, p):
        if da_fn(B, 1e6, *p) < 0.5:
            return np.nan
        try:
            return brentq(lambda F: da_fn(B, F, *p) - 0.5, 1e-6, 1e6, xtol=1e-8)
        except Exception:
            return np.nan
    def unpack(t):
        return tuple(np.exp(tv) if l else tv
                     for tv, l in zip(t, [n not in ("h1", "h2") for n in pnames]))
    def resid(t):
        p = unpack(t)
        ra = (da_fn(B_E, F_E, *p) - R_E) / np.std(R_E)
        rk = [np.log10(K12(B, p) / o) if np.isfinite(K12(B, p)) else 1.0
              for B, o in zip(K12_B, K12_obs)]
        return np.concatenate([ra, np.sqrt(lam) * np.array(rk) / 0.5])
    best = None
    for seed in range(nstart):
        r_ = np.random.default_rng(seed)
        xs = x0 + r_.normal(0, 0.7, len(x0))
        xs = np.clip(xs, np.asarray(lo) + 1e-9, np.asarray(hi) - 1e-9)
        try:
            sol = least_squares(resid, xs, bounds=(lo, hi), max_nfev=15000)
            ss = float(np.sum(resid(sol.x) ** 2))
            if best is None or ss < best[0]:
                best = (ss, sol.x)
        except Exception:
            pass
    if best is None:
        return float("nan"), float("nan"), tuple(float("nan") for _ in pnames)
    p = unpack(best[1])
    pr = da_fn(B_E, F_E, *p)
    r2a = 1 - ((pr - R_E) ** 2).sum() / ((R_E - R_E.mean()) ** 2).sum()
    errs = [abs(np.log10(K12(B, p) / o)) if np.isfinite(K12(B, p)) else 1.0
            for B, o in zip(K12_B, K12_obs)]
    return r2a, float(np.mean(errs)), p

def main():
    data = load_data()
    B_E, F_E, R_E, K12_B, K12_obs = data

    print("== §1 TCS depletion robustness ==")
    pn = ["Ki", "Ka", "N", "amax"]
    for rt, tag in [(RT, "R_T=17 nM"), (1e-9, "zero-depletion control")]:
        def da_rt(B, F, Ki, Ka, N, amax, _rt=rt):
            return da_dep(B, F, Ki, Ka, N, amax, _rt)
        r2, er, p = fit_joint(da_rt, pn, np.log([2, 200, 6, 1]),
                              np.log([0.05, 5, 1, 0.3]), np.log([300, 1e6, 60, 5]),
                              data, 1.0)
        print(f"  {tag}: λ=1 → R²={r2:.3f}, err={er:.3f} | "
              f"Ki={p[0]:.2f}, Ka={p[1]:.0f}, N={p[2]:.1f}, amax={p[3]:.2f}")

    print("\n== §2 gain retuning ==")
    pn2 = ["Ki0", "Ka", "N0", "Bg", "h1", "h2", "amax"]
    x0 = np.array([np.log(2), np.log(1e5), np.log(6), np.log(10), 0.5, 0.0, np.log(1)])
    lo = (np.log(0.05), np.log(5), np.log(1), np.log(0.01), -2.3, -2.3, np.log(0.3))
    hi = (np.log(300), np.log(1e6), np.log(60), np.log(1e4), 2.3, 2.3, np.log(5))
    for tag, fix in [("G-Ki", (None, 0.0)), ("G-N", (0.0, None)), ("G-both", (None, None))]:
        for lam in [1.0, 10.0, 30.0]:
            def da_v(B, F, Ki0, Ka, N0, Bg, h1, h2, amax,
                     _f=fix):
                h1 = _f[0] if _f[0] is not None else h1
                h2 = _f[1] if _f[1] is not None else h2
                return da_gain(B, F, Ki0, Ka, N0, Bg, h1, h2, amax)
            # with h fixed, still run full-parameter (h initialized at fixed value, bounds clamped)
            lo_v = list(lo); hi_v = list(hi)
            if fix[0] is not None: lo_v[4], hi_v[4] = fix[0] - 1e-6, fix[0] + 1e-6
            if fix[1] is not None: lo_v[5], hi_v[5] = fix[1] - 1e-6, fix[1] + 1e-6
            r2, er, p = fit_joint(da_v, pn2, x0, lo_v, hi_v, data, lam, nstart=8)
            print(f"  {tag} λ={lam:5.1f}: R²={r2:.3f}, err={er:.3f} | "
                  f"Ki0={p[0]:.2f}, Ka={p[1]:.0f}, N0={p[2]:.1f}, Bg={p[3]:.2f}, "
                  f"h1={p[4]:.2f}, h2={p[5]:.2f}, amax={p[6]:.2f}")

    print("""
== Summary ==
Depletion: front unmoved digit for digit → formally excluded.
Gain retuning: best-in-field K1/2 end (0.052 dex) but amplitude end collapses to 0.68,
  and N0, h hit non-physical bounds (gain rising without limit with background) → minimal form cannot close the fracture.
The fracture is robust to depletion / two subpopulations / incomplete adaptation / readout nonlinearity / gain retuning.
""")

if __name__ == "__main__":
    main()
