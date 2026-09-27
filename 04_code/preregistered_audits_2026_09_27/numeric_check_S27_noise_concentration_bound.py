# -*- coding: utf-8 -*-
"""
数值核查 S2.7 · 噪声集中界（工作笔记 §9 的代数防错，非注册性质）
=====================================================================
定理 N 三条的数值钉：
  N1 普适集中：失真率 D_n 围绕期望集中，P(|D-E|>=t) <= 2exp(-n t^2/2)
     （McDiarmid 有界差分，只需噪声逐坐标独立，不需要 margin 假设）
  N2 期望 = margin 泛函：E[D] = mean_pairs Phi(-|dx|/(sqrt(2)sigma))
     <= inf_z { 0.5*mu(delta<sqrt(2)sigma z) + Phi(-z) }
  N3 原子墙：margin 分布 mu 在 0 处有质量 rho 的原子（平台/精确打结）
     => 打结对贡献地板 ~rho^2（约定介于 rho^2/2 与 rho^2），任何方法不可救
     （这些对里本来就没有顺序信息——噪声恰好只杀死本来无序之处）
三场景：
  1) 平滑随机信号（无平台）：理论期望 vs 经验均值 + 集中度
  2) 40% 平台信号：非打结对吻合、打结对翻转率恒 1、总失真分解逐项吻合
  3) p53 工况锚点：sigma=15% 幅度 => 脉冲-基线翻转率 Phi(-4.71) ~ 1.2e-6
运行：python 数值核查_S27_噪声集中界.py   （秒级，纯 CPU）
"""
import numpy as np
from scipy.stats import norm

rng = np.random.default_rng(7)
TRIALS = 400
SIGMA = 0.05

def pair_data(x):
    n = len(x)
    iu = np.triu_indices(n, 1)
    st = np.sign(x[:, None] - x[None, :])[iu]
    d = np.abs(x[:, None] - x[None, :])[iu]
    return st, d

def empirical(x, sigma, trials=TRIALS):
    n = len(x)
    st, _ = pair_data(x)
    iu = np.triu_indices(n, 1)
    rates = np.empty(trials)
    for r in range(trials):
        y = x + rng.normal(0, sigma, n)
        sy = np.sign(y[:, None] - y[None, :])[iu]
        rates[r] = np.mean(sy != st)
    return rates

print("=" * 64)
print("场景 1：平滑随机信号（margin 分布无原子）")
print("=" * 64)
x1 = rng.uniform(0, 1, 300)
st1, d1 = pair_data(x1)
e1 = norm.cdf(-d1 / (np.sqrt(2) * SIGMA)).mean()
r1 = empirical(x1, SIGMA)
print(f"理论期望 {e1:.6f} | 经验均值 {r1.mean():.6f} | 经验std {r1.std():.2e}")
t = 0.02
print(f"N1 集中界 t=0.02: 上界 {2*np.exp(-300*t*t/2):.4f} | 经验越界率 {np.mean(np.abs(r1-r1.mean())>t):.4f}")

print()
print("=" * 64)
print("场景 2：40% 精确平台（margin 分布在 0 处有原子）")
print("=" * 64)
x2 = np.clip(rng.uniform(0, 1, 300), 0.4, 1.0)
st2, d2 = pair_data(x2)
tie = st2 == 0
rho2 = tie.mean()
n = len(x2)
iu = np.triu_indices(n, 1)
r_nt, r_t = [], []
for _ in range(TRIALS):
    y = x2 + rng.normal(0, SIGMA, n)
    sy = np.sign(y[:, None] - y[None, :])[iu]
    r_nt.append(np.mean(sy[~tie] != st2[~tie]))
    r_t.append(np.mean(sy[tie] != st2[tie]))
e_nt = norm.cdf(-d2[~tie] / (np.sqrt(2) * SIGMA)).mean()
print(f"打结对占比 {rho2:.4f}")
print(f"非打结对：理论 {e_nt:.6f} | 经验 {np.mean(r_nt):.6f}")
print(f"打结对：经验翻转率 {np.mean(r_t):.4f}（恒≈1，N3 地板）")
print(f"总失真分解 {(1-rho2)*e_nt + rho2:.4f} | 经验总计 {np.mean([ (1-rho2)*a + rho2*b for a,b in zip(r_nt,r_t) ]):.4f}")

print()
print("=" * 64)
print("场景 3：p53 工况锚点（sigma = 15% 幅度）")
print("=" * 64)
print(f"脉冲-基线翻转率 Phi(-1/(0.15*sqrt2)) = {norm.cdf(-1/(0.15*np.sqrt(2))):.2e}")
print("=> 实验级噪声下，跨量级比较基本不可摧毁；脆弱区只有精确打结对")
