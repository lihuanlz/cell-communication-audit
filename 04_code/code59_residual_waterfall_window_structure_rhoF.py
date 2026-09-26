#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code 59 residual waterfall + within-stimulus-window temporal structure + ρ(F) recheck v1.0.0
================================================
Executable corrected version of the three "zero-cost sprint" routes (protocol facts: 10 s stimulus window,
30 s gaps not imaged, no B=3/30 levels, measurement before adaptation — no "60 s residual" to speak of).

§1 amplitude residual waterfall (unit table + per-cell, two conventions, M2c λ=1 and standard amplitude-optimal models)
  Results: residuals small and centered at every level under both conventions (per-cell median −0.03…+0.01),
  no B=10 cliff, no mechanism-switch threshold. Only systematic structure: at the highest-F end of each row the model
  slightly overestimates (+0.04…+0.065; B=0 row F>2 bin per-cell median −0.097),
  i.e. data saturate earlier than the model. Per-cell spread (IQR≈0.15–0.2) far exceeds model error,
  and the |>0.15| fraction is actually lowest at B=10 (16.5%) — spread comes from cell heterogeneity,
  not model failure.
  ⚠️ dual-record: during development this code once fed total concentration sv to the model as foreground F (same-ligand
  files' s field includes background), causing spurious large residuals at B≥10; direct single-file single-level check
  (220302_FOV1, F=2: measured median 0.311 vs table 0.299) located it, then fixed to
  F=sv−B. Lesson: pooled pipeline results must be spot-checked at single points.

§2 K1/2 waterfall (the fracture's true location)
  model−data (dex, M2c λ=1): B=0:+0.295, 0.01:+0.142, 0.1:+0.271,
  0.3:+0.201, 1:+0.094, 10:−0.268, 100:+0.002.
  S-shaped mismatch: low-B plateau model ~4.0 vs measured 2.0–2.9 (overestimate), B=10 model
  7.5 vs measured 13.9 (underestimate), B=100 exact. The fracture lies entirely in the
  transition-zone shape of the K1/2 statistic (B∈[1,10]); the amplitude table itself is sound.

§3 within-stimulus-window temporal structure (honest version of the "adaptation-lag failure" test)
  per (B,F), linear slope of a over the last 6 stimulus columns (3 s):
  small positive slope at low F (+0.01/column, response not yet at plateau), decaying to ~0 as F grows;
  structure identical at B=0/10/100, no high-background-specific sustained rebound.
  → no B-dependent signal of adaptation intrusion within the 10 s window; "adaptation-lag failure" hypothesis not supported.

§4 ρ(F) recheck (premise test for route two)
  G10-scale pointwise inversion: ρ(F)=59(F=10)→16(F=80), seemingly decreasing with F.
  But direct fit with pure-competition MWC (effective ligand=F+ρ·L, L=10 µM L-Asp):
  constant ρ suffices — ρ≈24 (M2c anchor)/22 (amplitude-optimal anchor)/9 (Emonet anchor),
  residuals [−0.045,+0.006,+0.014,0,−0.002], only the F=5 point off by −0.045.
  → the "2× overestimate" is a scale artifact of G10 phenomenological-kernel inversion, not allosteric evidence;
    constant-ρ pure competition already fits; no ρ(F) correction term needed.
    (Conclusive implication: per-µM potency of L-Asp on Tar ≈ 22–24× meAsp,
     Ki_LAsp≈Ki_MeAsp/23 — literature-checkable.)

Data: Moore 2024 Dryad doi:10.5061/dryad.nvx0k6dzz (CC0)
run: python3 代码59_残差瀑布_窗内时间结构_与rhoF复核.py
dependencies: numpy, scipy, pandas; 代码54 in same directory (FILES/extraction logic reused)
"""

import importlib.util
import os
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.optimize import least_squares

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
    return np.log((1 + L / Ki) / (1 + L / Ka))

def m2c(B, F, Ki=14.53, Ka=1e6, N=11.8, beta=0.995, amax=1.63):
    lam_eff = LAM + (1 - beta) * N * gL(B, Ki, Ka)
    return amax * (1 / (1 + np.exp(lam_eff))
                   - 1 / (1 + np.exp(lam_eff + N * (gL(B + F, Ki, Ka) - gL(B, Ki, Ka)))))

def std(B, F, Ki=30.7, Ka=133.0, N=39.0, amax=1.39):
    return amax * (ASTAR - 1 / (1 + np.exp(LAM + N * (gL(B + F, Ki, Ka) - gL(B, Ki, Ka)))))

def iter_cells(bg_targets):
    for name, bg in FILES.items():
        B = 100.0 if name == "230831_FOV2" else float(bg)
        if B not in bg_targets:
            continue
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
            yield B, a, s

def cell_curves(bg_targets):
    out = defaultdict(list)
    for B, a, s in iter_cells(bg_targets):
        lev = defaultdict(list)
        for row in range(35):
            sr = s[row]
            sv = float(sr.max())
            stim = np.where(sr >= sv - 1e-9)[0]
            pre = np.arange(0, stim[0]) if len(stim) else None
            if pre is None or len(pre) < 4 or len(stim) < 6:
                continue
            F = sv - B                       # key: s includes background, foreground F=sv−B
            if F <= 0:
                continue
            lev[round(F, 4)].append(float(np.median(a[row, pre[-4:]]))
                                    - float(np.median(a[row, stim[-6:]])))
        good = {k: float(np.median(v)) for k, v in lev.items() if len(v) >= 4}
        if len(good) >= 3:
            for k, v in good.items():
                out[B].append((k, v))
    return out

def stim_slopes(bg_targets, fmin=0.4):
    out = defaultdict(list)
    for B, a, s in iter_cells(bg_targets):
        for row in range(35):
            sr = s[row]
            sv = float(sr.max())
            stim = np.where(sr >= sv - 1e-9)[0]
            F = sv - B
            if F <= fmin or len(stim) < 6:
                continue
            y = a[row, stim[-6:]]
            out[(B, round(F, 1))].append(float(np.polyfit(np.arange(6), y, 1)[0]))
    return out

def main():
    df = pd.read_csv(c54.CSV)
    B_E = df["B_uM"].to_numpy(float)
    F_E = df["F_uM"].to_numpy(float)
    R_E = df["R_a"].to_numpy(float)

    print("== §1a unit-table amplitude residuals (model−data) ==")
    for B in sorted(set(B_E)):
        m = B_E == B
        r2 = np.median(np.abs(m2c(B_E[m], F_E[m]) - R_E[m]))
        rs = np.median(np.abs(std(B_E[m], F_E[m]) - R_E[m]))
        print(f"  B={B:7.2f}: |M2c residual| median={r2:.3f}, |standard residual| median={rs:.3f}")

    print("\n== §1b per-cell amplitude residuals (data−M2c) ==")
    cur = cell_curves(set(B_E))
    for B in sorted(cur):
        pts = np.array(cur[B])
        res = pts[:, 1] - m2c(B, pts[:, 0])
        print(f"  B={B:7.2f}: n={len(res):4d}, median={np.median(res):+.3f}, "
              f"IQR=[{np.percentile(res, 25):+.3f},{np.percentile(res, 75):+.3f}], "
              f"|>0.15|={np.mean(np.abs(res) > 0.15):.1%}")
        for Flo, Fhi in [(0, 0.5), (0.5, 2), (2, 10), (10, 50)]:
            m = (pts[:, 0] > Flo) & (pts[:, 0] <= Fhi)
            if m.sum() > 5:
                print(f"      F∈({Flo},{Fhi}]: median={np.median(res[m]):+.3f} (n={m.sum()})")

    print("\n== §2 K1/2 waterfall (model−data, dex) ==")
    k12, cens, tot = c54.extract_K12()
    pred = {0.0: 4.01, 0.01: 4.02, 0.1: 4.05, 0.3: 4.11, 1.0: 4.34, 10.0: 7.51, 100.0: 120.16}
    for B in sorted(k12):
        if len(k12[B]) < 10:
            continue
        d = np.log10(pred[B] / np.array(k12[B]))
        print(f"  B={B:7.2f}: n={len(d):3d} (censored {cens.get(B, 0)}), median Δ={np.median(d):+.3f} dex")

    print("\n== §3 within-stimulus-window slope (a/0.5s) ==")
    sl = stim_slopes({0.0, 10.0, 100.0})
    for (B, F) in sorted(sl):
        v = np.array(sl[(B, F)])
        print(f"  B={B:6.1f} F={F:5.1f}: n={len(v):4d}, slope median={np.median(v):+.4f}, "
              f"IQR=[{np.percentile(v, 25):+.4f},{np.percentile(v, 75):+.4f}]")

    print("\n== §4 ρ(F) recheck: pure-competition MWC (L-Asp arm, five points F=5–80) ==")
    ROOTD = "/mnt/agents/output/03_细胞线3/公开数据/Moore2024_Chemotaxis_FRET"
    lev_all = defaultdict(list)
    for name in ["210721_FOV1", "210721_FOV2", "210818_FOV1", "210818_FOV2",
                 "210907_FOV1", "210908_FOV1", "210908_FOV2"]:
        d = sio.loadmat(os.path.join(ROOTD, name + ".mat"))["reorgData"]["resp_data"][0, 0]
        for ci in range(d.shape[1]):
            try:
                a = d["a"][0, ci].astype(float)
                s = d["s"][0, ci].astype(float)
            except Exception:
                continue
            if a.shape != (35, 20):
                continue
            lev = defaultdict(list)
            for row in range(35):
                sr = s[row]
                sv = float(sr.max())          # cross-ligand files: s records foreground only (meAsp)
                stim = np.where(sr >= sv - 1e-9)[0]
                pre = np.arange(0, stim[0]) if len(stim) else None
                if pre is None or len(pre) < 4 or len(stim) < 6:
                    continue
                lev[round(sv, 4)].append(float(np.median(a[row, pre[-4:]]))
                                         - float(np.median(a[row, stim[-6:]])))
            for k, v in lev.items():
                if len(v) >= 4:
                    lev_all[k].append(float(np.median(v)))
    Fs5 = np.array([5, 10, 20, 40, 80.0])
    obs5 = np.array([np.median(lev_all[f]) for f in Fs5])

    def da_comp(F, L, rho, Ki, Ka, N, amax):
        return amax * (ASTAR - 1 / (1 + np.exp(LAM + N * (gL(F + rho * L, Ki, Ka) - gL(rho * L, Ki, Ka)))))

    for Ki_a, Ka_a, N_a, tag in [(14.53, 1e6, 11.8, "M2c anchor"), (0.81, 200, 6, "Emonet anchor"),
                                 (30.7, 133, 39.0, "amplitude-optimal anchor")]:
        def rc(t):
            rho, amax = np.exp(t)
            return da_comp(Fs5, 10, rho, Ki_a, Ka_a, N_a, amax) - obs5
        best = None
        for seed in range(6):
            r_ = np.random.default_rng(seed)
            xs = np.log([20, 1.6]) + r_.normal(0, 0.6, 2)
            try:
                sol = least_squares(rc, xs, bounds=(np.log([0.1, 0.3]), np.log([300, 3])))
                ss = float(np.sum(rc(sol.x) ** 2))
                if best is None or ss < best[0]:
                    best = (ss, sol.x)
            except Exception:
                pass
        rho, amax = np.exp(best[1])
        pred = da_comp(Fs5, 10, rho, Ki_a, Ka_a, N_a, amax)
        print(f"  {tag}: ρ={rho:.1f}, amax={amax:.2f}, residuals={np.round(pred - obs5, 3)}")

    print("""
== Summary ==
Amplitude sound (no B=10 cliff); fracture lies entirely in the K1/2 transition-zone S-shaped mismatch (§2).
Within-window temporal structure shows no B dependence → adaptation-lag failure hypothesis not supported.
L-Asp arm: constant-ρ pure competition already fits (ρ≈22–24); the "2× overestimate" is a G10 scale artifact.
""")

if __name__ == "__main__":
    main()
