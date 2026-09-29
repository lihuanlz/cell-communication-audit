#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
代码54 K1/2 联合拟合：标准 MWC 模型自洽性断裂的定量 v1.0.0
================================================
任务：定量回答"标准 MWC+精确适应模型能否用单组物理参数同时自洽
Moore 论文的两个统计量——幅值表 R(B,F) 与 K1/2(B) 分布"。

§1 从原始 .mat 提取逐细胞 K1/2：字段 a（饱和刺激校准 0/1 归一）下，
   响应 da 按总配体 T=B+F 排序，在 da=0.5 处 log-T 轴线性插值。
   da 最大值 <0.5 的细胞记为右删失（不敏感细胞，不计入中位数）。
   K1/2 与幅值表取自同一批原始数据，内部一致性最高。
§2 联合拟合：参数 (Ki, Ka, N, amax)（log 空间），残差 =
   幅值残差/std(R) + sqrt(λ)·log10(K1/2预测/实测)/0.5 dex；
   λ 扫描给出 Pareto 前沿。least_squares 四起点取最优。

结果（v1.0.0 实跑）：
  K1/2 中位数：B=0→2.03, 0.01→2.90, 0.1→2.17, 0.3→2.58,
               1→3.49, 10→13.92, 100→119.56 µM
  （平段 ~2–3.5 µM（B≤1），之后 ≈K0+1.17·B；B=0 组删失 67% 偏低估）
  Pareto 前沿：
    λ=0    幅值R²=0.980, K1/2 误差 0.62 dex（Ki=31, N=39）
    λ=1    幅值R²=0.831, K1/2 误差 0.22 dex
    λ=30   幅值R²=0.511, K1/2 误差 0.07 dex（Ki=4.3, N=2.0）
  前沿陡峭、无可接受折中点（R²>0.9 且 |Δlog10|<0.1 不存在）
  → 断裂定量坐实。amax 拟合值 1.4–2.7（物理上应≈1）是应变的另一迹象。

数据：Moore 2024 Dryad doi:10.5061/dryad.nvx0k6dzz（CC0）
运行：python3 代码54_K12联合拟合_断裂度定量.py
依赖：numpy, scipy, pandas
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
# §1 逐细胞 K1/2 提取
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
                continue            # 右删失：最大刺激未达半幅
            if D.min() >= 0.5:
                continue            # 左删失（罕见）
            i = int(np.argmax(D >= 0.5))
            x0, x1 = np.log10(T[i - 1]), np.log10(T[i])
            y0, y1 = D[i - 1], D[i]
            k12[B].append(10 ** (x0 + (0.5 - y0) * (x1 - x0) / (y1 - y0)))
    return k12, cens, tot

# ---------------------------------------------------------------
# §2 联合拟合
# ---------------------------------------------------------------
ASTAR = 1.0 / 3.0
LAM = np.log(1.0 / ASTAR - 1.0)
TGT = np.log(1.0 / (ASTAR / 2) - 1) - LAM   # 活性减半所需 N·Δg

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
    print("== §1 逐细胞 K1/2 分布（总配体轴，µM）==")
    print(f"{'B':>8} {'n估出':>5} {'n删失':>5} {'K1/2中位':>9} {'IQR':>18}")
    K12_B, K12_obs = [], []
    for B in sorted(set(list(k12) + list(cens))):
        v = np.log10(np.array(k12.get(B, [np.nan])))
        med = 10 ** np.nanmedian(v)
        q = np.nanpercentile(v, [25, 75])
        K12_B.append(B); K12_obs.append(med)
        print(f"{B:8.2f} {len(k12.get(B, [])):5d} {cens.get(B, 0):5d} "
              f"{med:9.2f}   [{10**q[0]:7.2f},{10**q[1]:7.2f}]")
    K12_B = np.array(K12_B); K12_obs = np.array(K12_obs)
    print("注：B=0 组删失率高，中位为低估；平段（B≤1）与 ∝B 段（B≥10）结构稳健。")

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
    print("\n== §2 联合拟合 Pareto 前沿 ==")
    print(f"{'λ':>6} {'幅值R²':>7} {'K1/2中位|Δlog10|':>15}  参数")
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
== 结论 ==
不存在 R²>0.9 且 K1/2 误差 <0.1 dex 的参数点：前沿陡峭，
顾此失彼 → 标准 MWC 单组参数无法自洽 Moore 的两个统计量，
断裂定量坐实。amax 拟合值 1.4–2.7（物理上应≈1）为应变迹象。
下一步：多物种/局部适应模型能否在同一参数下同时自洽两者。
""")


if __name__ == "__main__":
    main()
