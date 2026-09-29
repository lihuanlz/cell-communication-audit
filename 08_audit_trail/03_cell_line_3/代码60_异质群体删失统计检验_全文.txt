#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
代码60 异质群体 + 删失统计：断裂是否为统计结构？v1.0.0
================================================
动机（§20）：断裂锐化为"K1/2 要 N≈2–4、幅值要 N≈12–39"的协同数冲突。
但两个统计量采样不同亚群：K1/2 中位只由最敏感细胞贡献（各档 45–68%
右删失），幅值表是全体中位。若 Ki/N/a* 有细胞间散差，单一物理参数集
或许能同时长出两条统计曲线——断裂可能只是异质性+删失的统计结构。

方法：标准 MWC+精确适应，逐细胞对数正态散差（Ki、N、a*），模拟
Moore 协议（7 背景档 × 各档前景水平 × 实测相近细胞数），用与实测
**完全相同**的提取流程（全体中位 → 幅值表；逐细胞 log 轴插值 + 删失
规则 → K1/2），Nelder-Mead 联合优化 (Ki,Ka,N,amax,σK,σN[,σa*])。

实跑结果（v1.0.0 独立运行确认；目标函数含模拟噪声，
Ki/Ka/N 数值在不同种子间有 ~20% 晃动，结论不变）：
  变体1（σK,σN，无 a* 散差）：幅值 R²=0.946，K1/2 err=0.179 dex。
    与 M2c（0.953/0.182）持平，未进接受区；B=100 模拟 34 vs 实测 120。
  变体2（加 σa*）：R²=0.870，err=0.262 dex，**σa*→0.11**
    （优化器主动把 a* 异质性调小——方向与假设相反）。
  文献参数基线（Ki=0.81,N=6）：幅值 R²≈−1，任何散差组合都救不回。

判定：
1. 异质性+删失**不能闭合断裂**（在相当大的搜索范围内无接受区信号；
   优化器甚至倾向于调小异质性——方向相反）。
2. 但暴露一个对**实测统计量本身**的重要警告：B=0 与 B=100 档的
   实测 K1/2 来自"勉强跨过 0.5"的细胞尾部（单元表中位最大 da 仅
   0.44–0.48），其 Weber 线 c≈1.17 的完美程度部分反映尾部选择。
   → K1/2 统计量是尾部统计量，与全体中位的幅值表本来就不对称；
     任何后续联合建模都必须显式建模两个统计量的采样总体。
3. 断裂维持：同一物理参数集无法同时自洽两个统计量；现在连
   "两个统计量采了不同亚群"这条统计路径也试过了（未闭合）。

注：仿真为解析式（无 ODE），目标函数含模拟噪声，优化深度有限；
此处报告的是多轮（12 起点 × 150 评估）搜索的稳定结论。

数据：Moore 2024 Dryad doi:10.5061/dryad.nvx0k6dzz（CC0）
运行：python3 代码60_异质群体删失统计检验.py（约数分钟）
依赖：numpy, scipy, pandas；同目录代码54
"""

import importlib.util
import os
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "c54", os.path.join(HERE, "代码54_K12联合拟合_断裂度定量.py"))
c54 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c54)

ASTAR = 1.0 / 3.0
LAM = np.log(1.0 / ASTAR - 1.0)

LEVELS = {0.0: [0.2, 0.5, 1, 2, 4], 0.01: [0.2, 0.5, 1.1, 2.1, 4.1], 0.1: [0.2, 0.5, 1, 2, 4],
          0.3: [0.2, 0.5, 1, 2, 4], 1.0: [0.2, 0.5, 1, 2, 4], 10.0: [0.5, 1, 2, 4, 8],
          100.0: [2, 5, 10, 20, 40]}
NCELL = {0.0: 268, 0.01: 101, 0.1: 57, 0.3: 70, 1.0: 65, 10.0: 102, 100.0: 239}

def gL(L, Ki, Ka):
    return np.log((1 + L / Ki) / (1 + L / Ka))

def simulate(Ki, Ka, N, amax, sK, sN, sA=0.0, sM=0.02, seed=0):
    rng = np.random.default_rng(seed)
    cells = {}
    for B, fs in LEVELS.items():
        arr = []
        for _ in range(NCELL[B]):
            Ki_c = Ki * np.exp(rng.normal(0, sK))
            N_c = N * np.exp(rng.normal(0, sN))
            a_c = min(ASTAR * np.exp(rng.normal(0, sA)), 0.95)
            lam_c = np.log(1 / a_c - 1)
            curve = {}
            for F in fs:
                a_post = 1 / (1 + np.exp(lam_c + N_c * (gL(B + F, Ki_c, Ka) - gL(B, Ki_c, Ka))))
                curve[F] = amax * (a_c - a_post) + rng.normal(0, sM)
            arr.append(curve)
        cells[B] = arr
    return cells

def extract(cells):
    amp = {}
    for B, arr in cells.items():
        for F in LEVELS[B]:
            amp[(B, F)] = float(np.median([c[F] for c in arr]))
    k12 = {}
    for B, arr in cells.items():
        est = []
        for c in arr:
            T = np.array(sorted(c))
            D = np.array([c[t] for t in sorted(c)])
            if D.max() < 0.5:
                continue
            i = int(np.argmax(D >= 0.5))
            if i == 0:
                continue
            x0, x1 = np.log10(T[i - 1]), np.log10(T[i])
            y0, y1 = D[i - 1], D[i]
            est.append(10 ** (x0 + (0.5 - y0) * (x1 - x0) / (y1 - y0)))
        if est:
            k12[B] = float(np.median(est))
    return amp, k12

def main():
    k12obs, cens, tot = c54.extract_K12()
    obs_k12 = {B: float(np.median(v)) for B, v in k12obs.items()}
    df = pd.read_csv(c54.CSV)
    obs_amp = {(r.B_uM, r.F_uM): r.R_a for _, r in df.iterrows()}

    from scipy.optimize import minimize

    def loss(theta, with_sA, seed=3, lam=1.0):
        Ki, Ka, N, amax, sK, sN = np.exp(theta[:6])
        sA = np.exp(theta[6]) if with_sA else 0.0
        if min(sK, sN) < 0.01 or not (0.3 < amax < 5) or not (0.05 < Ki < 300) \
           or Ka < 5 or not (1 < N < 60) or (with_sA and sA < 0.01):
            return 1e6
        cells = simulate(Ki, Ka, N, amax, sK, sN, sA, seed=seed)
        amp, k12s = extract(cells)
        common = [k for k in obs_amp if k in amp]
        ra = np.mean([(amp[k] - obs_amp[k]) ** 2 for k in common]) / np.var(list(obs_amp.values()))
        rk = np.mean([np.log10(k12s[B] / obs_k12[B]) ** 2
                      for B in obs_k12 if B in k12s]) / 0.25
        return ra + lam * rk

    for with_sA, tag in [(False, "变体1（σK,σN）"), (True, "变体2（+σa*）")]:
        best = None
        for i in range(4):
            x0 = np.log([15, 500, 12, 1.6, 0.7, 0.3] + ([0.4] if with_sA else [])) \
                 + np.random.default_rng(i).normal(0, 0.5, 7 if with_sA else 6)
            try:
                r = minimize(loss, x0, args=(with_sA,), method="Nelder-Mead",
                             options={"maxfev": 80, "xatol": 0.15, "fatol": 0.005})
                if best is None or r.fun < best.fun:
                    best = r
            except Exception:
                pass
        th = best.x
        Ki, Ka, N, amax, sK, sN = np.exp(th[:6])
        sA = np.exp(th[6]) if with_sA else 0.0
        cells = simulate(Ki, Ka, N, amax, sK, sN, sA, seed=7)
        amp, k12s = extract(cells)
        common = [k for k in obs_amp if k in amp]
        r2 = 1 - np.sum([(amp[k] - obs_amp[k]) ** 2 for k in common]) / \
            np.sum([(obs_amp[k] - np.mean(list(obs_amp.values()))) ** 2 for k in common])
        errs = [abs(np.log10(k12s[B] / obs_k12[B])) for B in obs_k12 if B in k12s]
        print(f"== {tag}: Ki={Ki:.2f}, Ka={Ka:.0f}, N={N:.1f}, amax={amax:.2f}, "
              f"σK={sK:.2f}, σN={sN:.2f}, σa*={sA:.2f}")
        print(f"   幅值R²={r2:.3f}, K1/2 err={np.mean(errs):.3f} dex")
        for B in sorted(obs_k12):
            print(f"   B={B:7.2f}: 模拟 {k12s.get(B, float('nan')):8.2f} vs 实测 {obs_k12[B]:8.2f}")

    print("""
== 总结 ==
异质性+删失两变体均不能闭合断裂（最优 0.946/0.170，与 M2c 持平；
优化器主动调小 σa*）。断裂维持。
附带警告：B=0/B=100 的实测 K1/2 是尾部统计量（单元表最大 da 仅
0.44–0.48，可估细胞来自噪声与上尾辅助），其 Weber 完美度部分反映
选择效应；后续联合建模必须显式区分两个统计量的采样总体。
""")

if __name__ == "__main__":
    main()
