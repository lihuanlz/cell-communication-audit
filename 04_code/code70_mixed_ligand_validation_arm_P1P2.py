#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code 70: mixed-ligand validation arm P-1/P-2 adjudication v1.0.0
================================================
Pre-registration file: 04_细胞线4/结果/预注册_混合配体验证臂_P1P2.md
(sealed on 2026-08-16 before data unblinding; this code strictly follows its pipeline and criteria)

Validation arm (Moore 2024 Dryad, downloaded this round after passing Anubis verification via browser):
  V1: 230822_FOV1, 230823_FOV1/2   background 10µM Lasp+100µM meAsp, foreground meAsp
  V2: 231017_FOV1/2, 231018_FOV1/2 background 100µM meAsp+0.1µM Lasp, foreground Lasp
Effective-background conversion (§18.3/§19.4, ρ=23): V1 B_eff=330, V2 B_eff=102.3 µM meAsp equivalents.

Predictions (fixed in the pre-registration):
  P-1a structure: R² in logF coordinates − R² in logr coordinates ≥ +0.3
  P-1b numeric: V1 da(F)=0.0677·ln(1+F/0.570); V2 da(F_L)=0.1682·ln(1+40.35·F_L)
        tolerance |observed median − prediction| ≤ 0.08 (all shared F levels)
  P-2  K1/2 (total effective axis): V1=388 µM (±0.15 dex; high censoring predicted if the F ceiling is insufficient);
        V2=122 µM equivalents = foreground 0.86 µM Lasp (±0.15 dex, estimable requires n≥20)
Pipeline identical to code 58 (a field, da = median of 4 pre-stimulus points − median of last 6 stimulus points;
per-cell K1/2: linear interpolation across 0.5 on the log10 axis, censored if D.max()<0.5).
Run: python3 代码70_混合配体验证臂_P1P2判决.py
================================================
"""
import os, warnings
from collections import defaultdict
import numpy as np
import scipy.io as sio

warnings.filterwarnings("ignore")
ROOT = "/mnt/agents/output/03_细胞线3/公开数据/Moore2024_Chemotaxis_FRET"

ARMS = {
    # s_sub: the amount to subtract from the s-field level to obtain the foreground F.
    # V1 foreground meAsp is the same ligand as the main background component → s records total meAsp (105=100+5), F=sv−100;
    # V2 foreground Lasp is a different ligand → s records the foreground (§18 data fact), but the background contains 0.1 Lasp,
    #   so a "total Lasp" ambiguity exists → compute both readings and compare goodness of fit (see §1b ambiguity check in the output).
    "V1": {"files": ["230822_FOV1", "230823_FOV1", "230823_FOV2"],
           "bg_txt": "10uM Lasp + 100uM meAsp", "B_eff": 330.0, "rho_fg": 1.0, "s_sub": 100.0},
    "V2": {"files": ["231017_FOV1", "231017_FOV2", "231018_FOV1", "231018_FOV2"],
           "bg_txt": "100uM meAsp + 0.1uM Lasp", "B_eff": 102.3, "rho_fg": 23.0, "s_sub": 0.0},
    "V2_alt": {"files": ["231017_FOV1", "231017_FOV2", "231018_FOV1", "231018_FOV2"],
           "bg_txt": "100uM meAsp + 0.1uM Lasp (s read as total Lasp)", "B_eff": 102.3, "rho_fg": 23.0, "s_sub": 0.1},
}
A_, Kc_, Bs_, p_ = 0.2268, 0.570, 69.7, 0.778   # G10 (P7 refit, pre-registered values)
RHO = 23.0
TOL_AMP, TOL_K12 = 0.08, 0.15

def g10(F_eff, B_eff):
    return A_ * np.log(1 + F_eff / Kc_) / (1 + (B_eff / Bs_) ** p_)

def unit_table(files, s_sub=0.0):
    lev_all, percell = defaultdict(list), []
    for name in files:
        p = os.path.join(ROOT, name + ".mat")
        d = sio.loadmat(p)["reorgData"]["resp_data"][0, 0]
        for ci in range(d.shape[1]):
            try:
                a = d["a"][0, ci].astype(float); s = d["s"][0, ci].astype(float)
            except Exception:
                continue
            if a.shape != (35, 20) or s.shape != (35, 20):
                continue
            lev = defaultdict(list)
            for row in range(35):
                sr = s[row]; sv = float(sr.max())
                stim = np.where(sr >= sv - 1e-9)[0]
                pre = np.arange(0, stim[0]) if len(stim) else None
                if pre is None or len(pre) < 4 or len(stim) < 6:
                    continue
                da = float(np.median(a[row, pre[-4:]])) - float(np.median(a[row, stim[-6:]]))
                lev[round(sv - s_sub, 4)].append(da)
            good = {k: float(np.median(v)) for k, v in lev.items() if len(v) >= 4}
            if len(good) >= 3:
                percell.append(good)
                for k, v in good.items():
                    lev_all[k].append(v)
    return lev_all, percell

def percell_k12(percell, B_eff, rho_fg):
    """K1/2 on the total effective axis; returns (array of estimable values, total cell count)"""
    out = []
    for good in percell:
        ks = sorted(good)
        T = B_eff + rho_fg * np.array(ks)
        D = np.array([good[k] for k in ks])
        if D.max() < 0.5:
            continue
        i = int(np.argmax(D >= 0.5))
        if i == 0:
            continue
        x0, x1 = np.log10(T[i - 1]), np.log10(T[i])
        y0, y1 = D[i - 1], D[i]
        out.append(10 ** (x0 + (0.5 - y0) * (x1 - x0) / (y1 - y0)))
    return np.array(out), len(percell)

print("== 0. File identity check (s-field levels vs BackgroundList conditions) ==")
for arm, cfg in ARMS.items():
    for name in cfg["files"]:
        d = sio.loadmat(os.path.join(ROOT, name + ".mat"))["reorgData"]["resp_data"][0, 0]
        s0 = d["s"][0, 0].astype(float)
        smm = d["s_min_max"][0, 0].astype(float)
        print(f"  {name}: s row levels={sorted(set(np.round(s0.max(axis=1),3)))[:8]}... "
              f"s_min_max range=[{smm.min():.3g},{smm.max():.3g}]  (background: {cfg['bg_txt']})")

tables = {arm: unit_table(cfg["files"], cfg.get("s_sub", 0.0)) for arm, cfg in ARMS.items()}

print("\n== 1. Amplitude unit table (foreground axis) and P-1b zero-fit reconciliation ==")
verdict_p1b = {}
for arm, cfg in ARMS.items():
    lev_all, percell = tables[arm]
    tag = "(ambiguity control, not adjudicated)" if arm.endswith("_alt") else ""
    print(f"\n[{arm}] {len(percell)} cells; B_eff={cfg['B_eff']} µM equivalents {tag}")
    devs = []
    for k in sorted(lev_all):
        obs = float(np.median(lev_all[k]))
        pred = g10(cfg["rho_fg"] * k, cfg["B_eff"])
        devs.append(abs(obs - pred))
        print(f"   F={k:8.3f}: observed da={obs:+.3f} (n={len(lev_all[k])})  predicted={pred:+.3f}  deviation={obs-pred:+.3f}")
    mx = max(devs)
    verdict_p1b[arm] = mx <= TOL_AMP
    print(f"   → max deviation {mx:.3f}, tolerance {TOL_AMP}: {'hit' if verdict_p1b[arm] else 'miss'}")

print("\n== 2. P-1a structural adjudication (coordinate test) ==")
for arm, cfg in ARMS.items():
    lev_all, _ = tables[arm]
    ks = sorted(lev_all)
    F = np.array(ks); R = np.array([np.median(lev_all[k]) for k in ks])
    F_eff = cfg["rho_fg"] * F
    def r2(xv, yv):
        A = np.vstack([xv, np.ones_like(xv)]).T
        c = np.linalg.lstsq(A, yv, rcond=None)[0]
        r = yv - A @ c
        return 1 - np.sum(r**2) / np.sum((yv - yv.mean())**2)
    r2_logF = r2(np.log10(F_eff), R)
    r2_logr = r2(np.log10((cfg["B_eff"] + F_eff) / cfg["B_eff"]), R)
    ok = (r2_logF - r2_logr) >= 0.3
    print(f"[{arm}] R²(logF)={r2_logF:.3f}  R²(logr)={r2_logr:.3f}  diff={r2_logF-r2_logr:+.3f}  "
          f"{'hit (≥+0.3)' if ok else 'miss'}")

print("\n== 3. P-2 per-cell K1/2 (total effective axis) ==")
for arm, cfg in ARMS.items():
    lev_all, percell = tables[arm]
    v, ntot = percell_k12(percell, cfg["B_eff"], cfg["rho_fg"])
    fmax = max(lev_all)
    pred_total = 1.17 * (1.95 + cfg["B_eff"])
    print(f"\n[{arm}] estimable {len(v)}/{ntot} cells (censored {100*(1-len(v)/ntot):.0f}%);"
          f"foreground F ceiling={fmax} (observed axis)")
    print(f"   pre-registered prediction: total-axis K1/2={pred_total:.1f} µM equivalents ±{TOL_K12} dex")
    if len(v) >= 20:
        med = float(np.median(v))
        dd = abs(np.log10(med / pred_total))
        print(f"   observed median={med:.1f} µM equivalents, IQR=[{np.percentile(v,25):.1f},{np.percentile(v,75):.1f}]")
        print(f"   deviation={np.log10(med/pred_total):+.3f} dex → {'hit' if dd<=TOL_K12 else 'miss'}")
        if cfg["rho_fg"] != 1:
            fg = (med - cfg["B_eff"]) / cfg["rho_fg"]
            print(f"   foreground-axis readout: F*={fg:.3f} µM Lasp (pre-registered point estimate 0.86)")
    else:
        print(f"   estimable cells <20 → per pre-registered rule: cannot adjudicate (if highly censored, qualitatively consistent with the V1 censoring prediction:"
              f"{'yes' if len(v)/ntot<0.5 else 'no'})")

print("\n== 4. Incidental observation (fracture side, not adjudicated): does the amplitude coordinate structure persist under mixed backgrounds ==")
print("(see the R² difference in section 2; a formal fracture test is left to the joint fit — this code sets no adjudication)")
