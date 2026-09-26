#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
代码57 TCS 耗竭稳健性检查 + 增益重调（候选①）联合检验 v1.0.0
================================================
两部分：
§1 TCS 耗竭修正稳健性检查。此前所有受体层模型均为零耗竭近似
   （g(L) 直接用总配体）。TCS 统一方程 ξ=p/(1-p)+p/κ，κ=Kd/(nR_T)。
   估算：最密集情形 1e9 cells/mL × 1e4 受体/细胞 → R_T≈17 nM，
   κ≈Ki/0.017 ≈ 48 ≫ 1（零耗竭区）。本代码不显式近似，而是解
   自由配体隐式方程 L_T = L_free + R_T·L_free/(Kd+L_free)（Kd≈Ki），
   用自由配体代入 g(L)，重扫 Pareto 前沿。
   实跑：λ=1 时 (R², err) = (0.948, 0.257)，与零耗竭对照完全相同，
   参数逐位一致 → 耗竭项正式排除为断裂来源（与 κ≫1 估算一致）。

§2 候选① 增益重调（有效参数随背景/甲基化态变化），三种最小参数化：
   G-Ki：Ki_eff(B)=Ki0·(1+B/Bg)^h1
   G-N： N_eff(B)=N0·(1+B/Bg)^h2
   G-双：两者同时。
   判别不变：R²>0.9 且 K1/2 误差<0.1 dex。

实跑结果（v1.0.0，分块验证确认：§1 用本文件函数独立复跑逐位一致；
§2 各变体经内核多轮 8–16 起点复跑，G-Ki 两轮完全一致，G-N/G-双
取全部运行最优——多起点优化对起点数敏感，起点少会掉进只顾幅值的
局部最优（R²=0.980/err=1.00），掩盖折中盆地，已双录）：
  G-Ki：λ=1 → 0.949/0.254；λ=30 → 0.663/0.064
        （低 λ 时 h→−2.3 下界、Bg→大：退化为常数 Ki；λ=3 逐 B 显示
         低平台 ~3.7 vs 实测 2.0–2.9，B=10/100 各欠 0.22–0.26 dex）
  G-N ：λ=1 → 0.948/0.247；λ=30 → 0.678/0.052
        （N0→60、h→2.3 双双撞界：拟合想要无限制随背景上升的增益，
         物理上不合理；Ka→9–30 µM 远离文献锚点）
  G-双：与 G-N 重合（h1→0，Ki 退回常数）→ 自由度浪费。

全模型 Pareto 对比（λ=1 折中点 / K1/2 端点）：
  标准 MWC        0.831/0.224 ↔ 0.511/0.074
  M2c 不完全适应   0.953/0.182 ↔ 0.843/0.089   ← 当前最优
  FRET 幂律读数    0.950/0.257 ↔ 0.758/0.100
  增益重调 G-N     0.948/0.247 ↔ 0.678/0.052

判定：
1. 耗竭：正式排除（前沿逐位不动）。
2. 增益重调最小形式：不能闭合断裂。它能把 K1/2 端推到 0.052 dex
   （全场最佳），但代价是幅值 R²=0.68 且参数撞非物理边界
   （增益随背景无限上升）。没有任何模型进入接受区。
3. 断裂对以下扰动全部稳健：TCS 耗竭、两亚群、不完全适应、
   读数非线性（压缩/扩张）、平滑增益重调。
   剩余路径：P8 独立数据外推；或直接接受"Moore 数据内部存在
   标准机制家族无法自洽的统计结构"作为本分支的终审结论并双录。

数据：Moore 2024 Dryad doi:10.5061/dryad.nvx0k6dzz（CC0）
运行：python3 代码57_TCS耗竭稳健性与增益重调检验.py
依赖：numpy, scipy, pandas；同目录代码54（数据提取复用）
"""

import importlib.util
import os
import warnings

import numpy as np
from scipy.optimize import brentq, least_squares

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "c54", os.path.join(HERE, "代码54_K12联合拟合_断裂度定量.py"))
c54 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c54)

ASTAR = 1.0 / 3.0
LAM = np.log(1.0 / ASTAR - 1.0)
RT = 0.017  # µM，最密集情形

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

# ---------- §1 耗竭 ----------
def da_dep(B, F, Ki, Ka, N, amax, RT=RT):
    Bf, TF = Lfree(B, Ki, RT), Lfree(B + F, Ki, RT)
    return amax * (ASTAR - 1 / (1 + np.exp(LAM + N * (gL(TF, Ki, Ka) - gL(Bf, Ki, Ka)))))

# ---------- §2 增益重调 ----------
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

    print("== §1 TCS 耗竭稳健性 ==")
    pn = ["Ki", "Ka", "N", "amax"]
    for rt, tag in [(RT, "R_T=17 nM"), (1e-9, "零耗竭对照")]:
        def da_rt(B, F, Ki, Ka, N, amax, _rt=rt):
            return da_dep(B, F, Ki, Ka, N, amax, _rt)
        r2, er, p = fit_joint(da_rt, pn, np.log([2, 200, 6, 1]),
                              np.log([0.05, 5, 1, 0.3]), np.log([300, 1e6, 60, 5]),
                              data, 1.0)
        print(f"  {tag}: λ=1 → R²={r2:.3f}, err={er:.3f} | "
              f"Ki={p[0]:.2f}, Ka={p[1]:.0f}, N={p[2]:.1f}, amax={p[3]:.2f}")

    print("\n== §2 增益重调 ==")
    pn2 = ["Ki0", "Ka", "N0", "Bg", "h1", "h2", "amax"]
    x0 = np.array([np.log(2), np.log(1e5), np.log(6), np.log(10), 0.5, 0.0, np.log(1)])
    lo = (np.log(0.05), np.log(5), np.log(1), np.log(0.01), -2.3, -2.3, np.log(0.3))
    hi = (np.log(300), np.log(1e6), np.log(60), np.log(1e4), 2.3, 2.3, np.log(5))
    for tag, fix in [("G-Ki", (None, 0.0)), ("G-N", (0.0, None)), ("G-双", (None, None))]:
        for lam in [1.0, 10.0, 30.0]:
            def da_v(B, F, Ki0, Ka, N0, Bg, h1, h2, amax,
                     _f=fix):
                h1 = _f[0] if _f[0] is not None else h1
                h2 = _f[1] if _f[1] is not None else h2
                return da_gain(B, F, Ki0, Ka, N0, Bg, h1, h2, amax)
            # 固定 h 时仍按全参数跑（h 初始即固定值、上下界夹紧）
            lo_v = list(lo); hi_v = list(hi)
            if fix[0] is not None: lo_v[4], hi_v[4] = fix[0] - 1e-6, fix[0] + 1e-6
            if fix[1] is not None: lo_v[5], hi_v[5] = fix[1] - 1e-6, fix[1] + 1e-6
            r2, er, p = fit_joint(da_v, pn2, x0, lo_v, hi_v, data, lam, nstart=8)
            print(f"  {tag} λ={lam:5.1f}: R²={r2:.3f}, err={er:.3f} | "
                  f"Ki0={p[0]:.2f}, Ka={p[1]:.0f}, N0={p[2]:.1f}, Bg={p[3]:.2f}, "
                  f"h1={p[4]:.2f}, h2={p[5]:.2f}, amax={p[6]:.2f}")

    print("""
== 总结 ==
耗竭：前沿逐位不动 → 正式排除。
增益重调：K1/2 端全场最佳（0.052 dex）但幅值端崩到 0.68，
  且 N0、h 撞非物理边界（增益随背景无限上升）→ 最小形式不能闭合断裂。
断裂对耗竭/两亚群/不完全适应/读数非线性/增益重调全部稳健。
""")

if __name__ == "__main__":
    main()
