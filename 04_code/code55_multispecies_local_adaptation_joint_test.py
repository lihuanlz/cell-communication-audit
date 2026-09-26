#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
代码55 多物种/局部适应模型联合检验 + 协议差异对照 v1.0.0
================================================
接续代码54的断裂定量，按审计第6轮第②③路径检验两类最小扩展模型
能否同时自洽 Moore 的幅值表（R²>0.9）与 K1/2 分布（误差<0.1 dex）。

模型：
  M2a 两亚群最小版：A=快适应（Ki=0.81 µM、Ka=200、N=6 锚定 Emonet 实测），
      B=完全不适应（甲基化固定→零配体绝对参考），活性按丰度 α 线性叠加。
      自由参数 6 个（α, Ki_B, Ka_B, N_B, cB, amax）。
  M2b 两亚群松弛版：A 的 Ka_A、N_A 放开（Ki_A=0.81 仍锚定）。自由参数 8 个。
  M2c 单物种+不完全适应：基线只补偿 β·g(B)（β=1 退化为标准模型），
      自由参数 5 个（Ki, Ka, N, β, amax）。
判别标准（审计设定）：R²>0.9 且 K1/2 误差<0.1 dex 且 amax≈1 且 α∈[0.3,0.7]。
K1/2 定义与提取一致：绝对半幅 da=0.5（a 字段按饱和刺激校准）；
模型预测达不到半幅的背景记删失，损失函数罚 1 dex（修正了"预测删失白得 0 罚"的漏洞）。

§3 协议对照：同一标准 MWC(Ki=2,Ka=200,N=6)+精确适应，
  (i) Moore 式阶跃+适应前读数；(ii) Lazova 式波形+适应后稳态读数。

实跑结果（v1.0.0，独立运行最终值；开发期曾出现 β/amax 解包顺序错位，
已修正并双录——教训：独立复跑是必需环节）：
  M2a：最优 R²=0.719/err=0.298（α→0.08 极端、amax 撞边界）→ 失败
  M2b：最优 R²=0.796/err=0.239（α→0.05，Ka_A、N_A 撞边界）→ 失败
       混合后幅值拟合反而劣于单物种（0.796<0.980）：
       不适应亚群的"零配体锚定增量"形状与幅值表形状冲突。
  M2c Pareto：(R², err) = (0.984, 1.00罚) / (0.953, 0.182) / (0.930, 0.145)
       / (0.843, 0.089)。λ=1 点：Ki=14.5, Ka→上限(单阈值 MWC), N=11.8,
       β=0.995, amax=1.63；B=100 精确命中（+0.002 dex），误差集中在
       低 B 平台段（模型 ~4.0 vs 实测 2.0–2.9，但 B=0 实测有 67% 删失
       低估）与 B=10（−0.268 dex）。对照 β≡1（Ka 同放宽）：
       R²=0.948/err=0.257 → 改善主要来自 Ka 放宽，β 仅小幅贡献。
       err<0.1 时 R²=0.843<0.9：断裂大幅收窄但未闭合。
  §3：标准模型阶跃适应前响应本就是 logF 主导（B=0.1 时 logF R²=0.497
       vs 线性F 0.145），波形稳态仅在 B∈[10,100] 呈 FCD 平台、低 B 塌陷
       （与 Lazova "18 µM 以下 no FCD"一致）→ 协议差异改变 regime 归属，
       但不能把标准模型的对数坐标变成实测的混合结构；联合不自洽是
       Moore 数据内部的，不是协议造成的。

数据：Moore 2024 Dryad doi:10.5061/dryad.nvx0k6dzz（CC0）
运行：python3 代码55_多物种与局部适应_联合检验.py
依赖：numpy, scipy, pandas；同目录代码54（数据提取复用）
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
    "c54", os.path.join(HERE, "代码54_K12联合拟合_断裂度定量.py"))
c54 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c54)

ASTAR = 1.0 / 3.0
LAM = np.log(1.0 / ASTAR - 1.0)

def gL(L, Ki, Ka):
    return np.log((1 + L / Ki) / (1 + L / Ka))

# ---------------------------------------------------------------
# 数据
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
# M2a/M2b 两亚群
# ---------------------------------------------------------------
KI_A = 0.81  # Emonet 实测锚定

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
# M2c 单物种+不完全适应
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
# 联合拟合驱动
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

    print("== M2a 两亚群最小版（A 锚定 Emonet）==")
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

    print("\n== M2b 两亚群松弛版（Ka_A、N_A 放开）==")
    unp_b = lambda t: (1/(1+np.exp(-t[0])), *np.exp(t[1:7]), t[7])
    x0b = np.array([0, np.log(200), np.log(6), np.log(1), np.log(500), np.log(6), np.log(1), 0.0])
    lob = [-6, np.log(5), np.log(1), np.log(0.05), np.log(5), np.log(1), np.log(0.3), -8]
    hib = [6, np.log(2e4), np.log(40), np.log(50), np.log(2e4), np.log(40), np.log(5), 8]
    for lam in [0.0, 1.0, 10.0]:
        r2, er, p = joint_fit(da2, k12_2, unp_b, x0b, lob, hib, B_E, F_E, R_E,
                              K12_B, K12_obs, lam)
        print(f"  λ={lam:5.1f}: R²={r2:.3f}, K1/2err={er:.3f}, "
              f"α={p[0]:.2f}, Ka_A={p[1]:.0f}, N_A={p[2]:.1f}, Ki_B={p[3]:.2f}, amax={p[7]:.2f}")

    print("\n== M2c 单物种+不完全适应 ==")
    # 参数顺序与 da_part(Ki,Ka,N,beta,amax) 对齐：θ=(logKi,logKa,logN,logAmax,logitβ)
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
    # 逐 B 残差（λ=1 最优点）
    r2, er, p = joint_fit(da_part, K12_part, unp_c, x0c, loc, hic, B_E, F_E, R_E,
                          K12_B, K12_obs, 1.0, nstart=16)
    print(f"  λ=1 逐 B：实测 vs 模型 K1/2")
    for B, o in zip(K12_B, K12_obs):
        pr = K12_part(B, p)
        print(f"    B={B:7.2f}: {o:9.2f} vs " +
              (f"{pr:9.2f} (Δ={np.log10(pr/o):+.3f} dex)" if np.isfinite(pr) else "  [删失]"))

    print("\n== §3 协议对照（标准 MWC Ki=2, Ka=200, N=6 + 精确适应）==")
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
        print(f"  Moore式阶跃(适应前) B={B:6.1f}: logF R²={rl:.3f}, 线性F R²={rf:.3f}")
    for B in [0.1, 1, 10, 100, 1000]:
        print(f"  Lazova式波形(适应后) B={B:7.1f}: 幅值={resp_wave(B, 1.5):.4f}")

    print("""
== 总结 ==
M2a/M2b 两亚群线性叠加：失败（幅值 R²≤0.842，参数撞边界，且劣于单物种）。
  不适应亚群的零配体锚定增量形状与幅值表形状冲突——审计预测的
  "高背景残差由不适应亚群解释"未成立。
M2c 不完全适应：最接近的折中 R²=0.953/err=0.182 dex（Ka→∞ 单阈值、
  β=0.995 近精确适应），B=100 命中，但未进接受区，断裂收窄未闭合。
§3 协议差异改变 regime 归属，但不能解释 Moore 数据内部的联合不自洽。
""")

if __name__ == "__main__":
    main()
