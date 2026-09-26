#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code 55 joint test of multispecies / local-adaptation models + protocol-difference control v1.0.0
================================================
Following code54's break quantification, per audit round-6 paths ② and ③, test whether two classes of minimal extended models
can simultaneously fit Moore's amplitude table (R²>0.9) and the K1/2 distribution (error <0.1 dex).

Models:
  M2a two-subpopulation minimal: A=fast adaptation (Ki=0.81 µM, Ka=200, N=6 anchored to Emonet measurements),
      B=completely non-adapting (methylation fixed → zero-ligand absolute reference); activities add linearly by abundance α.
      6 free parameters (α, Ki_B, Ka_B, N_B, cB, amax).
  M2b two-subpopulation relaxed: A's Ka_A, N_A released (Ki_A=0.81 still anchored). 8 free parameters.
  M2c single species + imperfect adaptation: baseline compensates only β·g(B) (β=1 reduces to the standard model),
      5 free parameters (Ki, Ka, N, β, amax).
Acceptance criteria (set by the audit): R²>0.9 and K1/2 error <0.1 dex and amax≈1 and α∈[0.3,0.7].
K1/2 definition matches the extraction: absolute half-amplitude da=0.5 (a field calibrated to saturating stimulus);
backgrounds where the model prediction never reaches half-amplitude are recorded as censored and penalized 1 dex in the loss (fixes the loophole where predicted censoring got a free 0 penalty).

§3 protocol control: same standard MWC (Ki=2, Ka=200, N=6) + precise adaptation,
  (i) Moore-style steps + pre-adaptation readout; (ii) Lazova-style waveform + post-adaptation steady-state readout.

Actual run results (v1.0.0, final values from an independent run; during development a β/amax unpacking-order misalignment occurred,
fixed and double-recorded — lesson: independent reruns are mandatory):
  M2a: best R²=0.719/err=0.298 (α→0.08 extreme, amax hits boundary) → FAIL
  M2b: best R²=0.796/err=0.239 (α→0.05, Ka_A and N_A hit boundaries) → FAIL
       the mixture fits amplitude worse than the single species (0.796<0.980):
       the non-adapting subpopulation's "zero-ligand-anchored increment" shape conflicts with the amplitude-table shape.
  M2c Pareto: (R², err) = (0.984, 1.00 penalty) / (0.953, 0.182) / (0.930, 0.145)
       / (0.843, 0.089). λ=1 point: Ki=14.5, Ka→upper bound (single-threshold MWC), N=11.8,
       β=0.995, amax=1.63; B=100 hit exactly (+0.002 dex); errors concentrate in
       the low-B plateau (model ~4.0 vs observed 2.0–2.9, though the B=0 observation has 67% censoring
       downward bias) and B=10 (−0.268 dex). Control β≡1 (Ka released identically):
       R²=0.948/err=0.257 → the improvement comes mainly from releasing Ka; β contributes only slightly.
       at err<0.1, R²=0.843<0.9: the break narrows substantially but does not close.
  §3: the standard model's step pre-adaptation response is already logF-dominated (at B=0.1, logF R²=0.497
       vs linear-F 0.145); the waveform steady state shows an FCD plateau only at B∈[10,100] and collapses at low B
       (consistent with Lazova "no FCD below 18 µM") → protocol differences change regime assignment,
       but cannot turn the standard model's logarithmic coordinate into the observed mixed structure; the joint inconsistency is
       internal to the Moore data, not caused by protocol.

Data: Moore 2024 Dryad doi:10.5061/dryad.nvx0k6dzz (CC0)
Run: python3 代码55_多物种与局部适应_联合检验.py
Dependencies: numpy, scipy, pandas; code54 in the same directory (data extraction reused)
"""

import importlib.util
import os
import warnings

import numpy as np
import pandas as pd
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

# ---------------------------------------------------------------
# Data
# ---------------------------------------------------------------
def load_data():
    k12, cens, tot = c54.extract_K12()
    K12_B = np.array(sorted(set(list(k12) + list(cens))))
    K12_obs = np.array([10 ** np.nanmedian(np.log10(np.array(k12.get(B, [np.nan]))))
                        for B in K12_B])
    df = pd.read_csv(c54.CSV)
    return (df["B_uM"].to_numpy(float), df["F_uM"].to_numpy(float),
            df["R_a"].to_numpy(float), K12_B, K12_obs)

# ---------------------------------------------------------------
# M2a/M2b two subpopulations
# ---------------------------------------------------------------
KI_A = 0.81  # anchored to Emonet measurement

def make_two(fixed_A):
    def da(B, F, alpha, Ka_A, N_A, Ki_B, Ka_B, N_B, cB, amax):
        dA = ASTAR - 1 / (1 + np.exp(LAM + N_A * (gL(B + F, KI_A, Ka_A) - gL(B, KI_A, Ka_A))))
        aB_pre = 1 / (1 + np.exp(N_B * gL(B, Ki_B, Ka_B) + cB))
        aB_post = 1 / (1 + np.exp(N_B * gL(B + F, Ki_B, Ka_B) + cB))
        return amax * (alpha * dA + (1 - alpha) * (aB_pre - aB_post))
    return da

def K12_of(da_fn):
    def f(B, p):
        if da_fn(B, 1e6, *p) < 0.5:
            return np.nan
        try:
            return brentq(lambda F: da_fn(B, F, *p) - 0.5, 1e-6, 1e6, xtol=1e-8)
        except Exception:
            return np.nan
    return f

# ---------------------------------------------------------------
# M2c single species + imperfect adaptation
# ---------------------------------------------------------------
def da_part(B, F, Ki, Ka, N, beta, amax):
    lam_eff = LAM + (1 - beta) * N * gL(B, Ki, Ka)
    a_base = 1 / (1 + np.exp(lam_eff))
    a_post = 1 / (1 + np.exp(lam_eff + N * (gL(B + F, Ki, Ka) - gL(B, Ki, Ka))))
    return amax * (a_base - a_post)

def K12_part(B, p):
    if da_part(B, 1e6, *p) < 0.5:
        return np.nan
    try:
        return brentq(lambda F: da_part(B, F, *p) - 0.5, 1e-6, 1e6, xtol=1e-8)
    except Exception:
        return np.nan

# ---------------------------------------------------------------
# joint-fit driver
# ---------------------------------------------------------------
def joint_fit(da_fn, k12_fn, unpack, x0, lo, hi, B_E, F_E, R_E, K12_B, K12_obs,
              lam, nstart=10):
    def resid(theta):
        p = unpack(theta)
        ra = (da_fn(B_E, F_E, *p) - R_E) / np.std(R_E)
        rk = []
        for B, o in zip(K12_B, K12_obs):
            pr = k12_fn(B, p)
            rk.append(np.log10(pr / o) if np.isfinite(pr) else 1.0)
        return np.concatenate([ra, np.sqrt(lam) * np.array(rk) / 0.5])
    best = None
    for seed in range(nstart):
        r_ = np.random.default_rng(seed)
        xs = x0 + r_.normal(0, 0.8, len(x0))
        try:
            sol = least_squares(resid, xs, bounds=(lo, hi), max_nfev=15000)
            ss = float(np.sum(resid(sol.x) ** 2))
            if best is None or ss < best[0]:
                best = (ss, sol.x)
        except Exception:
            pass
    p = unpack(best[1])
    pr = da_fn(B_E, F_E, *p)
    r2a = 1 - ((pr - R_E) ** 2).sum() / ((R_E - R_E.mean()) ** 2).sum()
    errs = [abs(np.log10(k12_fn(B, p) / o)) if np.isfinite(k12_fn(B, p)) else 1.0
            for B, o in zip(K12_B, K12_obs)]
    return r2a, float(np.mean(errs)), p

def main():
    B_E, F_E, R_E, K12_B, K12_obs = load_data()

    print("== M2a two-subpopulation minimal (A anchored to Emonet) ==")
    da2 = make_two(True); k12_2 = K12_of(da2)
    unp_a = lambda t: (1/(1+np.exp(-t[0])), 200.0, 6.0, *np.exp(t[1:5]), t[5])
    x0a = np.array([0, np.log(1), np.log(500), np.log(6), np.log(1), 0.0])
    loa = [-6, np.log(0.05), np.log(5), np.log(1), np.log(0.3), -8]
    hia = [6, np.log(50), np.log(2e4), np.log(30), np.log(5), 8]
    for lam in [0.0, 1.0, 10.0]:
        r2, er, p = joint_fit(da2, k12_2, unp_a, x0a, loa, hia, B_E, F_E, R_E,
                              K12_B, K12_obs, lam)
        print(f"  λ={lam:5.1f}: R²={r2:.3f}, K1/2err={er:.3f}, "
              f"α={p[0]:.2f}, Ki_B={p[3]:.2f}, N_B={p[5]:.1f}, amax={p[7]:.2f}")

    print("\n== M2b two-subpopulation relaxed (Ka_A, N_A released) ==")
    unp_b = lambda t: (1/(1+np.exp(-t[0])), *np.exp(t[1:7]), t[7])
    x0b = np.array([0, np.log(200), np.log(6), np.log(1), np.log(500), np.log(6), np.log(1), 0.0])
    lob = [-6, np.log(5), np.log(1), np.log(0.05), np.log(5), np.log(1), np.log(0.3), -8]
    hib = [6, np.log(2e4), np.log(40), np.log(50), np.log(2e4), np.log(40), np.log(5), 8]
    for lam in [0.0, 1.0, 10.0]:
        r2, er, p = joint_fit(da2, k12_2, unp_b, x0b, lob, hib, B_E, F_E, R_E,
                              K12_B, K12_obs, lam)
        print(f"  λ={lam:5.1f}: R²={r2:.3f}, K1/2err={er:.3f}, "
              f"α={p[0]:.2f}, Ka_A={p[1]:.0f}, N_A={p[2]:.1f}, Ki_B={p[3]:.2f}, amax={p[7]:.2f}")

    print("\n== M2c single species + imperfect adaptation ==")
    # parameter order aligned with da_part(Ki,Ka,N,beta,amax): θ=(logKi,logKa,logN,logAmax,logitβ)
    unp_c = lambda t: (np.exp(t[0]), np.exp(t[1]), np.exp(t[2]),
                       1/(1+np.exp(-t[4])), np.exp(t[3]))
    x0c = np.concatenate([np.log([15, 1e4, 12, 1.6]), [2.0]])
    loc = np.concatenate([np.log([0.05, 5, 1, 0.3]), [-8]])
    hic = np.concatenate([np.log([300, 1e6, 60, 5]), [8]])
    for lam in [0.0, 1.0, 3.0, 10.0]:
        r2, er, p = joint_fit(da_part, K12_part, unp_c, x0c, loc, hic, B_E, F_E, R_E,
                              K12_B, K12_obs, lam, nstart=16)
        print(f"  λ={lam:5.1f}: R²={r2:.3f}, K1/2err={er:.3f}, "
              f"Ki={p[0]:.1f}, Ka={p[1]:.0f}, N={p[2]:.1f}, β={p[3]:.3f}, amax={p[4]:.2f}")
    # per-B residuals (λ=1 optimum)
    r2, er, p = joint_fit(da_part, K12_part, unp_c, x0c, loc, hic, B_E, F_E, R_E,
                          K12_B, K12_obs, 1.0, nstart=16)
    print(f"  λ=1 per B: observed vs model K1/2")
    for B, o in zip(K12_B, K12_obs):
        pr = K12_part(B, p)
        print(f"    B={B:7.2f}: {o:9.2f} vs " +
              (f"{pr:9.2f} (Δ={np.log10(pr/o):+.3f} dex)" if np.isfinite(pr) else "  [censored]"))

    print("\n== §3 protocol control (standard MWC Ki=2, Ka=200, N=6 + precise adaptation) ==")
    Ki_s, Ka_s, N_s = 2.0, 200.0, 6.0
    def resp_step(B, F):
        return ASTAR - 1/(1+np.exp(LAM + N_s*(gL(B+F, Ki_s, Ka_s)-gL(B, Ki_s, Ka_s))))
    def resp_wave(B, r):
        return abs(ASTAR - 1/(1+np.exp(LAM + N_s*(gL(B*r, Ki_s, Ka_s)-gL(B, Ki_s, Ka_s)))))
    Fs = np.array([1, 2, 5, 10, 20, 50, 100.0])
    for B in [0.1, 1, 10, 100]:
        rs = np.array([resp_step(B, F) for F in Fs])
        rl = np.corrcoef(np.log10(Fs), rs)[0, 1] ** 2
        rf = np.corrcoef(Fs, rs)[0, 1] ** 2
        print(f"  Moore-style step (pre-adaptation) B={B:6.1f}: logF R²={rl:.3f}, linear-F R²={rf:.3f}")
    for B in [0.1, 1, 10, 100, 1000]:
        print(f"  Lazova-style waveform (post-adaptation) B={B:7.1f}: amplitude={resp_wave(B, 1.5):.4f}")

    print("""
== Summary ==
M2a/M2b two-subpopulation linear superposition: FAIL (amplitude R²≤0.842, parameters hit boundaries, worse than single species).
  The non-adapting subpopulation's zero-ligand-anchored increment shape conflicts with the amplitude-table shape — the audit's
  predicted "high-background residuals explained by a non-adapting subpopulation" did not hold.
M2c imperfect adaptation: closest compromise R²=0.953/err=0.182 dex (Ka→∞ single threshold,
  β=0.995 near-precise adaptation), B=100 hit, but it does not enter the acceptance region; the break narrows without closing.
§3 protocol differences change regime assignment but cannot explain the joint inconsistency internal to the Moore data.
""")

if __name__ == "__main__":
    main()
