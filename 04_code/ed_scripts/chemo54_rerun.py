"""Re-run code 54 logic on archived Moore 2024 .mat data:
per-cell K1/2 extraction + joint-fit Pareto sweep. Saves JSON for EDFig. 10."""
import os, json, warnings
from collections import defaultdict
import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.optimize import brentq, least_squares
warnings.filterwarnings("ignore")

WS = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911"
ROOT = os.path.join(WS, r"02_细胞线2\公开数据\Moore2024_Chemotaxis_FRET")
CSV = os.path.join(WS, r"03_细胞线3\结果\P7_Moore2024_FCD-Weber\代码49b_单元表.csv")
OUT = os.path.join(WS, r"05_主线纲领与设计\论文_细胞通讯审计_2026-09-23\figures_svg\ED\scripts\chemo54_rerun.json")

FILES = {
    "210802_FOV1": 0, "210802_FOV2": 0, "210805_FOV1": 0, "210805_FOV2": 0,
    "220106_FOV1": 0, "230417_FOV1": 0,
    "230815_FOV1": 0.01, "230815_FOV2": 0.01, "230816_FOV1": 0.01, "230816_FOV2": 0.01,
    "230830_FOV1": 0.1, "230830_FOV2": 0.1, "230831_FOV1": 0.1, "230831_FOV2": 0.1,
    "220615_FOV1": 0.3, "230410_FOV1": 0.3, "230428_FOV1": 1.0, "230429_FOV1": 1.0,
    "220302_FOV1": 10.0, "220303_FOV1": 10.0,
    "210816_FOV1": 100.0, "210816_FOV2": 100.0, "230717_FOV1": 100.0, "230718_FOV1": 100.0,
}

def extract_K12():
    k12 = defaultdict(list); cens = defaultdict(int); tot = defaultdict(int)
    for name, bg in FILES.items():
        p = os.path.join(ROOT, name + ".mat")
        if not os.path.isfile(p):
            print("missing", name); continue
        try:
            rd = sio.loadmat(p)["reorgData"]["resp_data"][0, 0]
        except Exception as e:
            print("load fail", name, e); continue
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
                continue
            if D.min() >= 0.5:
                continue
            i = int(np.argmax(D >= 0.5))
            x0, x1 = np.log10(T[i - 1]), np.log10(T[i])
            y0, y1 = D[i - 1], D[i]
            k12[B].append(10 ** (x0 + (0.5 - y0) * (x1 - x0) / (y1 - y0)))
    return k12, cens, tot

ASTAR = 1.0 / 3.0
LAM = np.log(1.0 / ASTAR - 1.0)
TGT = np.log(1.0 / (ASTAR / 2) - 1) - LAM

def gL(L, Ki, Ka):
    return np.log((1 + L / Ki) / (1 + L / Ka))

def amp_model(B, F, Ki, Ka, N, amax):
    return amax * (ASTAR - 1 / (1 + np.exp(LAM + N * (gL(B + F, Ki, Ka) - gL(B, Ki, Ka)))))

def K12_model(B, Ki, Ka, N):
    try:
        return brentq(lambda F: N * (gL(B + F, Ki, Ka) - gL(B, Ki, Ka)) - TGT, 1e-6, 1e6, xtol=1e-8)
    except Exception:
        return np.nan

def main():
    k12, cens, tot = extract_K12()
    K12_B, K12_obs, K12_iqr, K12_all = [], [], [], {}
    print("== per-background K1/2 ==")
    for B in sorted(set(list(k12) + list(cens))):
        v = np.log10(np.array(k12.get(B, [np.nan])))
        med = 10 ** np.nanmedian(v)
        q = np.nanpercentile(v, [25, 75])
        K12_B.append(B); K12_obs.append(med); K12_iqr.append([10 ** q[0], 10 ** q[1]])
        K12_all[str(B)] = [float(x) for x in k12.get(B, [])]
        print(f"B={B:8.2f} n_est={len(k12.get(B, [])):4d} n_cens={cens.get(B,0):4d} "
              f"med={med:8.2f} IQR=[{10**q[0]:.2f},{10**q[1]:.2f}]")
    K12_B = np.array(K12_B); K12_obs = np.array(K12_obs)

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
    front = []
    print("== Pareto sweep ==")
    for lam in [0.0, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0]:
        best = None
        for seed in range(4):
            r_ = np.random.default_rng(seed)
            xs = x0 + r_.normal(0, 0.8, 4)
            try:
                sol = least_squares(resid, xs, args=(lam,),
                                    bounds=(np.log([0.05, 5, 2, 0.3]), np.log([300, 2e4, 60, 3])),
                                    max_nfev=8000)
                ss = np.sum(resid(sol.x, lam) ** 2)
                if best is None or ss < best[0]:
                    best = (ss, sol.x)
            except Exception:
                pass
        Ki, Ka, N, amax = np.exp(best[1])
        pr = amp_model(B_E, F_E, Ki, Ka, N, amax)
        r2a = 1 - ((pr - R_E) ** 2).sum() / ((R_E - R_E.mean()) ** 2).sum()
        dk = float(np.mean([abs(np.log10(K12_model(B, Ki, Ka, N) / o)) for B, o in zip(K12_B, K12_obs)]))
        print(f"lam={lam:5.1f} R2={r2a:.3f} dK={dk:.3f} Ki={Ki:.1f} Ka={Ka:.0f} N={N:.1f} amax={amax:.2f}")
        front.append(dict(lam=lam, r2=float(r2a), dk=dk, Ki=float(Ki), Ka=float(Ka), N=float(N), amax=float(amax)))

    json.dump(dict(K12_B=list(map(float, K12_B)), K12_obs=list(map(float, K12_obs)),
                   K12_iqr=K12_iqr, K12_all=K12_all,
                   cens={str(k): int(v) for k, v in cens.items()},
                   tot={str(k): int(v) for k, v in tot.items()},
                   front=front),
              open(OUT, "w"), indent=1)
    print("saved", OUT)

if __name__ == "__main__":
    main()
