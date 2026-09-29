#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
代码56 FRET 读数层非线性检验 v1.0.0
================================================
候选②：坐标变换是否发生在读数层。上游固定为标准 MWC+精确适应
（Ki, Ka, N 自由），CheY-P 线性于活性，检验三种 CheY-P→FRET 映射：
  V1 饱和结合映射 S=Y/(1+Y)，Y=c·a（CheY-P/CheZ 结合的双曲线），绝对响应；
  V2 同上映射 + 按各背景饱和校准脉冲归一（对应 a 字段的校准方式）；
  V3 幂律映射 S∝Y^q（q>1 扩张 / q<1 压缩），检验坐标变换方向。
判别同前：幅值 R²>0.9 且 K1/2 误差<0.1 dex。K1/2 为绝对半幅 da=0.5，
模型预测删失罚 1 dex。least_squares 14 起点取最优。

实跑结果（v1.0.0）：
  V1 饱和映射：c→0.34–0.37（CheY-P 始终处于读数近线性区，饱和不被
    数据需要），amax 撞边界 5。λ=1：R²=0.946/err=0.256。
    → 饱和方向（压缩）与需求相反：实测比标准模型"更增量"，
      压缩只会让 logF 更像 logF。
  V2 按背景归一：更差（λ=1：R²=0.754/err=0.206，c→0.05 撞下界）。
  V3 幂律映射：q→2.0–2.1（扩张）稳定出现，λ=1：R²=0.950/err=0.257，
    amax 撞边界 5。优于标准 MWC（0.831/0.224 在同样 λ 下幅值更高），
    但不优于 M2c 不完全适应+单阈值（0.953/0.182），且未进接受区。
判定：读数层非线性三种形态均不能闭合断裂。数据需要的方向是"扩张"
  （q≈2 被稳定选出），但单独读数层变换不够；M2c 仍是当前最优
  （R²=0.953/err=0.182，Ka→∞ 单阈值对数感受器）。
意义：候选②基本排除为唯一主因。断裂的闭合需要受体层本身具有
  "低 B 平台 K1/2≈2.5 µM + B=10 陡起到 14 µM"的结构——指向候选①
  （有效参数随甲基化态/背景变化）或读数扩张与受体对数化的组合。

数据：Moore 2024 Dryad doi:10.5061/dryad.nvx0k6dzz（CC0）
运行：python3 代码56_FRET读数层非线性检验.py
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

def da_satn(B, F, Ki, Ka, N, c, amax):         # V2（amax 占位不用）
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
    for name, fn, qmode in [("V1 饱和映射", da_sat, False),
                            ("V2 按背景归一", da_satn, False),
                            ("V3 幂律映射", da_pow, True)]:
        print(f"== {name} ==")
        lo_q = 0.2 if qmode else 0.05
        for lam in [0.0, 1.0, 10.0]:
            if qmode:
                # V3 参数序 (Ki,Ka,N,q,amax)
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
== 总结 ==
V1/V2 饱和映射：被数据拒绝（c→小、近线性区；压缩方向与需求相反）。
V3 幂律映射：q≈2（扩张）稳定出现，λ=1 达 R²=0.950/err=0.257，
  优于纯标准 MWC 但不优于 M2c（0.953/0.182），未进接受区。
读数层非线性作为唯一主因基本排除；数据需要的是"扩张"方向的变换，
断裂闭合指向受体层结构（候选①增益重调）或与读数扩张的组合。
""")

if __name__ == "__main__":
    main()
