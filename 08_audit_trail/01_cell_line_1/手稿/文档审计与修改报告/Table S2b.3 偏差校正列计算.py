# -*- coding: utf-8 -*-
"""
Table S2b.3 偏差校正列计算（审计补充脚本，2026-07-31）
=====================================================
用途：为 SI Table S2b.3 新增"校正后偏差"列，基于 SI S2b.6-7 optimality and
LoD tails.py 的同一套精确引擎（无 Monte Carlo），逐行计算：

  bias_exact   bias(exact)             R1 估计器的精确偏差（%）
  bias_pred    bias(1/2 f''Var pred.)  二阶 delta 法预测偏差（%）
  diff(oracle) = bias_exact - bias_pred 相减列（oracle：f'' 取在真值 P_pos 处）
  bias_corr    逐试验精确校正后的残余偏差（可实施口径，等价于原脚本 TABLE C：
               M_corr = M_hat - 0.5*f''(P_hat)*Var(P_hat)，f'' 逐试验取值）

重要结论（2026-07-31 复算）：
  A 列（相减列）与 B 列（逐试验精确校正）在 M 较大时不一致——
  f'' 是凸函数，E[f''(P_hat)] > f''(E[P_hat])，可实施校正减得更多。
  正文 "subtracting the term gives a bias-corrected estimator" 对应的是 B 列。
  建议表中新列用 B 列数值；masked 两行（M=500/800, kappa=1）填 n/a。

几何：panel 1（Omega=1000, N=n=50, beta=20, b=0.01），与原脚本一致。
已验证：kappa=1 时 M=50/100/200 的 bias_corr = -0.03/-0.15/-2.41，
与原脚本 TABLE C 输出完全一致。
运行时间约 1-2 分钟（TABLE A 引擎为主）。仅需 numpy/scipy。
"""

import numpy as np
from scipy.special import gammaln
from scipy.stats import binom

np.seterr(all='ignore')

# ============================================================
# 精确引擎（与 SI S2b.6-7 optimality and LoD tails.py 相同）
# ============================================================

def canonical_probs(W, Omega, KV):
    """P(C|W), canonical ensemble (S1c.3), log-shift normalised."""
    maxC = min(W, Omega); C = np.arange(maxC + 1)
    lw = (gammaln(W + 1) - gammaln(C + 1) - gammaln(W - C + 1)
          + gammaln(Omega + 1) - gammaln(C + 1) - gammaln(Omega - C + 1)
          + gammaln(C + 1) - C * np.log(KV))
    lw -= lw.max(); w = np.exp(lw)
    return C, w / w.sum()


def poisson_grid(M):
    """W grid and Poisson(M) weights, tail-safe (14 sigma)."""
    Wmax = int(M + 14 * np.sqrt(M + 1)) + 2
    W = np.arange(0, Wmax + 1)
    lw = W * np.log(M) - gammaln(W + 1) - M
    pw = np.exp(lw - lw.max())
    return W, pw / pw.sum()


def alloc_DP(Omega, N, beta, Cmax):
    """B[s, C] = [t^C]((1+t)^beta - 1)^s * 2^(-beta*s), exact zeros masked."""
    k = np.arange(0, min(beta, Cmax) + 1)
    kern = np.exp(gammaln(beta + 1) - gammaln(k + 1) - gammaln(beta - k + 1)) * 2.0 ** (-beta)
    kern[0] = 0.0
    B = np.zeros((N + 1, Cmax + 1)); B[0, 0] = 1.0
    for s in range(1, N + 1):
        B[s] = np.convolve(B[s - 1], kern)[:Cmax + 1]
    return B


def PZ_given_W(W, Omega, N, beta, KV, B, n_obs):
    """P(Z_sig = z | W), z = 0..n_obs, marginalised over C|W."""
    C, w = canonical_probs(W, Omega, KV)
    C = C.astype(int)
    logC_Om = gammaln(Omega + 1) - gammaln(C + 1) - gammaln(Omega - C + 1)
    PzC = np.zeros((len(C), n_obs + 1))
    logC_N = gammaln(N + 1) - gammaln(np.arange(N + 1) + 1) - gammaln(N - np.arange(N + 1) + 1)
    for s in range(0, N + 1):
        BsC = B[s, C]
        pos = BsC > 0
        if not np.any(pos):
            continue
        logQ = np.full(len(C), -np.inf)
        logQ[pos] = np.log(BsC[pos]) + beta * s * np.log(2) + logC_N[s] - logC_Om[pos]
        Qs = np.exp(np.clip(logQ, -745, 50))
        zmax = min(s, n_obs)
        zs = np.arange(0, zmax + 1)
        hw = np.exp(gammaln(s + 1) - gammaln(zs + 1) - gammaln(s - zs + 1)
                    + gammaln(N - s + 1) - gammaln(n_obs - zs + 1) - gammaln(N - s - n_obs + zs + 1)
                    - (gammaln(N + 1) - gammaln(n_obs + 1) - gammaln(N - n_obs + 1)))
        PzC[:, :zmax + 1] += Qs[:, None] * hw[None, :]
    return np.sum(w[:, None] * PzC, axis=0)


def PZ_and_score(M, Omega, N, beta, kappa, n_obs, b=0.0, B=None):
    """Exact P(Z=z) and dP(Z=z)/dM (Poisson score), z = 0..n_obs."""
    if B is None:
        B = alloc_DP(Omega, N, beta, min(Omega, int(M + 14 * np.sqrt(M + 1)) + 2))
    W, pw = poisson_grid(M); KV = kappa * Omega
    Pz = np.zeros(n_obs + 1); Sz = np.zeros(n_obs + 1)
    for i, wv in enumerate(W):
        if pw[i] < 1e-15:
            continue
        pzw = PZ_given_W(int(wv), Omega, N, beta, KV, B, n_obs)
        Pz += pw[i] * pzw
        Sz += pw[i] * (wv / M - 1.0) * pzw
    if b > 0:  # Z = Z_sig + Binomial(n_obs - Z_sig, b)
        PzB = np.zeros(n_obs + 1); SzB = np.zeros(n_obs + 1)
        for s in range(n_obs + 1):
            if Pz[s] < 1e-18 and abs(Sz[s]) < 1e-18:
                continue
            ks = np.arange(0, n_obs - s + 1)
            bm = binom.pmf(ks, n_obs - s, b)
            PzB[s:s + len(ks)] += Pz[s] * bm
            SzB[s:s + len(ks)] += Sz[s] * bm
        Pz, Sz = PzB, SzB
    return Pz / Pz.sum(), Sz / Pz.sum()


def fisher_exact(M, Omega, N, beta, kappa, n_obs, b=0.0, B=None):
    """I(Z;M) = sum_z (dP/dM)^2 / P; returns (I, Pz)."""
    Pz, Sz = PZ_and_score(M, Omega, N, beta, kappa, n_obs, b, B)
    return np.sum(Sz ** 2 / np.maximum(Pz, 1e-300)), Pz


def ppos_exact(M, Omega, N, beta, kappa, b=0.0):
    """P_pos(M) = b + (1-b)(1-eps_1), exact."""
    KV = kappa * Omega
    W, pw = poisson_grid(M)
    e1 = 0.0
    for i, wv in enumerate(W):
        if pw[i] < 1e-15:
            continue
        C, w = canonical_probs(int(wv), Omega, KV)
        r = np.exp(gammaln(Omega - beta + 1) - gammaln(Omega - beta - C + 1)
                   - gammaln(Omega + 1) + gammaln(Omega - C + 1))
        e1 += pw[i] * np.sum(w * r)
    return b + (1 - b) * (1 - e1)


def r1_invert(Ph, kappa, Omega, beta, b):
    """R1 inversion P_hat -> M_hat (master equation), with background correction."""
    Ps = (Ph - b) / (1 - b)
    if Ps <= 0 or Ps >= 1:
        return np.nan
    x = 1 - Ps
    pp = 1 - x ** (1.0 / beta)
    return (pp / (1 - pp) + pp / kappa) * kappa * Omega


# ============================================================
# Geometry: panel 1（与 Table S2b.3 一致）
# ============================================================
Om1, N1, beta1, n1, b1 = 1000, 50, 20, 50, 0.01

print("Building allocation DP (panel-1 geometry)...")
B1 = alloc_DP(Om1, N1, beta1, Om1)
print("done.\n")

print(f"{'M':>4s} {'kap':>4s} {'bias_exact%':>11s} {'bias_pred%':>10s} "
      f"{'diff(oracle)%':>13s} {'bias_corr_exact%':>16s}")

# Table S2b.3 的 9 行
rows = [(10, 0.1), (50, 0.1), (100, 0.1),
        (10, 1.0), (50, 1.0), (100, 1.0), (200, 1.0), (500, 1.0), (800, 1.0)]

h = 1e-4
for M, kap in rows:
    I, Pz = fisher_exact(M, Om1, N1, beta1, kap, n1, b=b1, B=B1)
    zs = np.arange(n1 + 1)
    Mh = np.array([r1_invert(z / n1, kap, Om1, beta1, b1) for z in zs])
    okm = np.isfinite(Mh) & (Pz > 1e-15)
    P_ok = Pz[okm] / Pz[okm].sum()

    # bias (exact)
    EM = np.sum(P_ok * Mh[okm])
    bias_ex = (EM - M) / M * 100

    # bias (1/2 f''Var pred.)：f'' 取在真值 P_pos 处
    Ppos = ppos_exact(M, Om1, N1, beta1, kap, b=b1)
    f2 = (r1_invert(Ppos + h, kap, Om1, beta1, b1)
          - 2 * r1_invert(Ppos, kap, Om1, beta1, b1)
          + r1_invert(Ppos - h, kap, Om1, beta1, b1)) / h ** 2
    Ph = zs[okm] / n1
    VPh = np.sum(P_ok * Ph ** 2) - np.sum(P_ok * Ph) ** 2
    bias_pred = 0.5 * f2 * VPh / M * 100

    # 逐试验精确校正（TABLE C 口径）：f'' 在每个 P_hat 上取值
    f2h = np.array([(r1_invert(p + h, kap, Om1, beta1, b1)
                     - 2 * r1_invert(p, kap, Om1, beta1, b1)
                     + r1_invert(p - h, kap, Om1, beta1, b1)) / h ** 2 for p in Ph])
    Mc = Mh[okm] - 0.5 * f2h * VPh
    bias_corr = (np.sum(P_ok * Mc) - M) / M * 100

    print(f"{M:4d} {kap:4.1f} {bias_ex:11.2f} {bias_pred:10.2f} "
          f"{bias_ex - bias_pred:13.2f} {bias_corr:16.2f}")

print()
print("说明：masked 两行 (M=500/800, kap=1) 在 Table S2b.3 中填 n/a；")
print("建议新列采用 bias_corr_exact（逐试验精确校正，可实施口径）。")
