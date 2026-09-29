# -*- coding: utf-8 -*-
"""
代码83：为什么数字化、为什么噪声在上游 —— L1 损失分解 + 编码器结构变分竞赛
=================================================================================
日期：2026-09-25 ｜ 种子固定：20260925 ｜ 确定性部分全程确定

谱系：
  - 编码器层封卷：判词卡/封卷卡 2026-09-18（代码51r3/51r4A 归档参数）
  - 全链路级联审计：代码81（判词卡 2026-09-23）——本脚本直接复用其全部常数与互信息口径
  - 机制模型审计：代码82（判词卡 2026-09-23）——"94.2/5.5 是噪声标定口径函数"的登记出处

问题（用户命题，2026-09-25）：
  Q-A：代码81 判得 L1 传感层损失 0.872 bits = 31% H(D) 是全链路最弱层。
       这 0.872 里多少是物理上不可约的（DSB Poisson 涨落 + 修复稀释），
       多少是映射设计选择（聚类常数 c=50）？最优量化器能挽回多少？
  Q-B：把同一份上游噪声（Poisson DSB + 聚类 + Exp(1.0T) 首脉冲弥散 + IPI 抖动 CV 0.085）
       平等注入三种候选编码器结构：
         臂A 固定周期计数（观察到的架构：T=5.5h 钉死、幅度钉死、剂量走计数）
         臂B 周期调制（FM：剂量走周期 T(D)，脉冲数固定）
         臂C 幅度调制（AM：剂量走幅度 A(D)，脉冲数固定）
       在各统计量"天然承受"的散布下（计数：上游注入；周期：归档 CV(T) 钉死；
       幅度：细胞间增益 LogNormal sigma_k），谁的端到端信息最大？
       若观察到的架构=竞赛赢家，则"细胞选数字化"从归纳升级为变分最优性陈述。

公平性约束（双录）：
  - 三臂共用同一 L0/L1 上游生成器与同口径 t1 弥散；
  - 臂B/C 脉冲数固定 N=5（隔离编码变量，防止计数信息泄漏进 FM/AM 臂）；
  - 臂B 的 T(D) 动态范围（2 倍）与臂C 的 A(D) 动态范围（约 2 倍）对齐；
  - 解码级（B2）三臂共用同一对泄漏积分器（p21 型 tau=10h / PUMA 型 tau=4h）。
"""

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path
import sys
from math import lgamma

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
try:
    from daimon_runtime import setup_plot
    setup_plot()
except Exception:
    pass

SEED = 20260925
ROOT = Path(__file__).resolve().parent.parent
OUT_JSON = ROOT / "结果" / "代码83_变分竞赛_结果.json"
OUT_PNG = ROOT / "结果" / "代码83_变分竞赛_六面体.png"
OUT_SVG = ROOT / "结果" / "代码83_变分竞赛_六面体.svg"
plt.rcParams["svg.fonttype"] = "none"

# ---------------- 与代码81逐字一致的归档常数 ----------------
DOSES = np.array([0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0])   # Gy，等权 -> H(D)=log2(7)
N_CELLS = 4000
K0_DSB = 35.0
C_CLUSTER = 50.0
N_MAX = 12
T_PERIOD = 5.5
MU_D = 1.0 * T_PERIOD
JIT_IPI = 0.085
A0 = 1.0
A_CV = 0.002
W = 48.0
Q21, TAU21 = 100.0, 10.0
QPU, TAUPU, PUMA_NTHR = 100.0, 4.0, 3
H_D = float(np.log2(len(DOSES)))
N_BOOT = 200

# 臂B/C 专用常数（显式声明，非拟合）
N_FIXED = 5            # FM/AM 臂固定脉冲数（预算匹配：臂A 全体平均脉冲数 ~2.5，FM/AM 用 5，属对 FM/AM 有利的偏置，双录）
GAMMA_T = 0.15         # T(D) = T_PERIOD * (D/1Gy)^(-GAMMA_T)：D 0.1->10 对应 T 7.77->3.89 h（2 倍范围）
AMP_SLOPE = 0.15       # A(D) = A0*(1+0.15*ln(D/1Gy))（代码81 假想幅度臂口径）
AMP_INTR = 0.15        # 假想幅度臂内生 CV（代码81 同口径）
ETA_OBS = 0.02

# ---------------- 互信息工具（与代码81同口径） ----------------

def mi_discrete(x, y):
    x = np.asarray(x); y = np.asarray(y)
    Kx = int(x.max()) + 1; Ky = int(y.max()) + 1
    N = len(x)
    cxy = np.bincount(x * Ky + y, minlength=Kx * Ky).reshape(Kx, Ky) / N
    cx = cxy.sum(1, keepdims=True); cy = cxy.sum(0, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        term = cxy * np.log2(cxy / (cx @ cy))
    return float(np.nansum(term))


def mi_mm(x, y):
    x = np.asarray(x); y = np.asarray(y)
    rx = len(np.unique(x)); ry = len(np.unique(y)); N = len(x)
    bias = (rx - 1) * (ry - 1) / (2 * N * np.log(2))
    return max(0.0, mi_discrete(x, y) - bias)


def qbin(v, nbin):
    v = np.asarray(v, dtype=float)
    edges = np.unique(np.quantile(v, np.linspace(0, 1, nbin + 1)))
    if len(edges) <= 2:
        return np.zeros(len(v), dtype=int)
    return np.clip(np.digitize(v, edges[1:-1]), 0, len(edges) - 2)

# ---------------- Part A：L1 损失分解 ----------------

def log_poisson_pmf(n, lam):
    n = np.asarray(n, dtype=float)
    return n * np.log(lam) - lam - np.vectorize(lgamma)(n + 1)


def poisson_channel(doses, k0, phi=1.0, vmax=600):
    """P(n|D) 矩阵，n=0..vmax；phi 为存活（未修复）比例。"""
    ns = np.arange(vmax + 1)
    P = np.zeros((len(doses), vmax + 1))
    for di, D in enumerate(doses):
        lam = k0 * D * phi
        P[di] = np.exp(log_poisson_pmf(ns, lam))
    return P / P.sum(1, keepdims=True), ns


def mi_of_quantizer(P, assign):
    """P: (nDose, nVal) 条件分布；assign: 每个取值 -> 箱号。返回 I(D; bin)。
    下溢防护：pbD < 1e-15 的贡献视为 0（其真实贡献 < 1e-13 bits，可忽略）。"""
    nD = P.shape[0]
    prior = 1.0 / nD
    nb = assign.max() + 1
    pbD = np.zeros((nb, nD))
    for b in range(nb):
        m = assign == b
        pbD[b] = P[:, m].sum(1) * prior
    p_b = pbD.sum(1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.clip(pbD / (p_b * prior), 1e-300, 1e300)
        t = np.where(pbD > 1e-15, pbD * np.log2(ratio), 0.0)
    return float(np.nansum(t))


def optimal_quantizer_dp(P, M):
    """精确 DP：连续箱划分最小化 sum_b p(b) H(D|b)（等价最大化 I(D;bin)）。
    返回 (assign, I_opt)。复杂度 O(V^2 M)，V<=601。"""
    nD, V = P.shape
    prior = 1.0 / nD
    # bin(i..j) 的成本：p(bin) * H(D|bin)，用累计和 O(1) 算
    cs = np.cumsum(P, axis=1)  # (nD, V)
    neg_inf = -1e18
    # 预计算所有 bin 成本
    cost = np.full((V + 1, V + 1), np.inf)  # cost[i][j] = bin i..j-1
    for i in range(V):
        pj = np.zeros(nD)
        cum_prev = cs[:, i - 1] if i > 0 else np.zeros(nD)
        # 逐 j 扩展太慢则向量化：直接矩阵算
    # 向量化版：对每对 (i,j)
    idx = np.arange(V)
    for i in range(V):
        seg = cs[:, i:] - (cs[:, i - 1:i] if i > 0 else np.zeros((nD, 1)))  # (nD, V-i)
        pbin = seg.sum(0) * prior  # (V-i,)
        with np.errstate(divide="ignore", invalid="ignore"):
            post = seg * prior / pbin[None, :]
            ent = -np.nansum(np.where(post > 0, post * np.log2(post), 0.0), axis=0)
        cost[i, i + 1:] = pbin * ent
    # DP
    dp = np.full((M + 1, V + 1), np.inf)
    bp = np.zeros((M + 1, V + 1), dtype=int)
    dp[0, 0] = 0.0
    for m in range(1, M + 1):
        for j in range(1, V + 1):
            # dp[m][j] = min_i dp[m-1][i] + cost[i][j]
            cand = dp[m - 1, :j] + cost[:j, j]
            k = int(np.argmin(cand))
            dp[m, j] = cand[k]
            bp[m, j] = k
    # 回溯
    bounds = []
    j = V
    for m in range(M, 0, -1):
        i = bp[m, j]
        bounds.append(i)
        j = i
    bounds = sorted(set(bounds))
    assign = np.zeros(V, dtype=int)
    b = 0
    for i in range(V):
        while b < len(bounds) and i >= bounds[b]:
            b += 1
        assign[i] = b
    I_opt = H_D - dp[M, V]
    return assign, float(I_opt), [int(x) for x in bounds]


def part_A():
    log = lambda *a: print(*a, flush=True)
    log("[A] L1 损失分解启动")
    P, ns = poisson_channel(DOSES, K0_DSB)
    I_DSB = mi_of_quantizer(P, ns)  # 恒等映射（每取值一箱）
    log(f"  I(D;N_DSB) 精确值 = {I_DSB:.4f} bits（代码81 MC 口径 2.526）")

    # A1：聚类常数扫描
    C_LIST = [1, 2, 5, 10, 25, 50, 100, 200]
    c_sweep = []
    for c in C_LIST:
        assign = np.minimum(np.ceil(ns / c).astype(int), N_MAX)
        # ceil(0/c)=0 -> 与 1..c 同箱？ceil(0)=0 单独一箱，保持
        I_c = mi_of_quantizer(P, assign)
        c_sweep.append(dict(c=c, I=I_c, loss=I_DSB - I_c))
        log(f"  c={c:>4}: I(D;N_trig)={I_c:.4f} bits, 映射损失={I_DSB - I_c:.4f}")

    # A2：M=13 级最优量化器（脉冲帽约束下的信息上限）
    assign_opt, I_opt, bounds = optimal_quantizer_dp(P, N_MAX + 1)
    log(f"  最优 {N_MAX + 1} 级量化器: I={I_opt:.4f} bits, 损失={I_DSB - I_opt:.4f}, "
        f"边界 N={bounds}")

    # A3：修复稀释（thinning）——传感层"物理上"要付的部分
    PHI = [1.0, 0.8, 0.6, 0.4, 0.25, 0.15, 0.1, 0.05]
    thin_sweep = []
    for phi in PHI:
        Pt, _ = poisson_channel(DOSES, K0_DSB, phi=phi)
        _, ns_t = poisson_channel(DOSES, K0_DSB, phi=phi)
        ns_t = np.arange(Pt.shape[1])
        I_full = mi_of_quantizer(Pt, ns_t)
        # 同一 phi 下的最优 13 级量化
        ao, Io, _ = optimal_quantizer_dp(Pt, N_MAX + 1)
        # 同一 phi 下的 c=50 聚类
        ac = np.minimum(np.ceil(ns_t / C_CLUSTER).astype(int), N_MAX)
        Ic = mi_of_quantizer(Pt, ac)
        thin_sweep.append(dict(phi=phi, I_perfect_sense=I_full,
                               I_optimal_quant=Io, I_c50=Ic,
                               loss_thinning=I_DSB - I_full,
                               loss_total_optimal=I_DSB - Io,
                               loss_total_c50=I_DSB - Ic))
        log(f"  phi={phi:.2f}: 完美传感 I={I_full:.4f} | 最优量化 {Io:.4f} | "
            f"c=50 聚类 {Ic:.4f} | 稀释损失 {I_DSB - I_full:.4f}")

    return dict(I_DSB_exact=I_DSB, L0_poisson_floor=H_D - I_DSB,
                c_sweep=c_sweep, optimal_quantizer=dict(M=N_MAX + 1, I=I_opt,
                loss=I_DSB - I_opt, bounds=bounds),
                thinning_sweep=thin_sweep)

# ---------------- Part B：编码器结构竞赛 ----------------

def gen_train(n_trig, T_code, n, rng):
    """脉冲串生成：t1~Exp(1.0*T_cell)，IPI 抖动 CV 0.085，T_cell 逐细胞散布由调用方给。
    n_trig: (n,) 每细胞脉冲数；T_code: (n,) 编码周期（含逐细胞散布后的由调用方先乘好）。"""
    times = np.full((n, N_MAX), np.nan)
    t1 = rng.exponential(T_code)  # 上游弥散口径：mu_d = 1.0 * T_cell（代码81: Exp(1.0*T)）
    tcur = t1.copy()
    for j in range(N_MAX):
        times[:, j] = tcur
        interval = T_code * (1.0 + JIT_IPI * rng.standard_normal(n))
        interval = np.maximum(interval, 0.3 * T_code)
        tcur = tcur + interval
    in_win = times <= W
    pidx = np.arange(N_MAX)[None, :]
    active = (pidx < n_trig[:, None]) & in_win
    return times, active


def decoder(times, active):
    """p21/PUMA 泄漏积分器（代码81 同口径）。"""
    with np.errstate(invalid="ignore"):
        m_p21 = np.nansum(np.where(active, Q21 * TAU21 * (1.0 - np.exp(-(W - times) / TAU21)), 0.0), axis=1)
        pidx = np.arange(N_MAX)[None, :]
        pm = active & (pidx >= PUMA_NTHR - 1)
        m_puma = np.nansum(np.where(pm, QPU * TAUPU * (1.0 - np.exp(-(W - times) / TAUPU)), 0.0), axis=1)
    return m_p21, m_puma


def part_B():
    log = lambda *a: print(*a, flush=True)
    log("[B] 编码器结构竞赛启动")
    rng = np.random.default_rng(SEED + 10)
    n_tot = N_CELLS * len(DOSES)
    d = np.repeat(np.arange(len(DOSES)), N_CELLS)
    D = DOSES[d]

    # 共同上游（三臂同一生成器同一随机流）
    n_dsb = rng.poisson(K0_DSB * D, n_tot)
    n_trig = np.minimum(np.ceil(n_dsb / C_CLUSTER).astype(int), N_MAX)

    # ---- 臂A：固定周期计数（观察到的架构） ----
    timesA, actA = gen_train(n_trig, np.full(n_tot, T_PERIOD), n_tot, rng)
    n_pulse_A = actA.sum(1)
    I_A_enc = mi_discrete(d, n_pulse_A)
    m21A, mpuA = decoder(timesA, actA)
    I_A_dec = mi_discrete(d, qbin(m21A, 12) * 12 + qbin(mpuA, 12))
    log(f"  臂A 计数: I(D;N)={I_A_enc:.4f} bits, 过解码器 I={I_A_dec:.4f}, "
        f"平均脉冲 {n_pulse_A.mean():.2f}")

    # ---- 臂B：周期调制（FM） ----
    T_code = T_PERIOD * np.power(D, -GAMMA_T)  # D=1 -> 5.5h
    CVT_LIST = [0.0, 0.016, 0.05, 0.1, 0.2]
    B_rows = []
    for cvT in CVT_LIST:
        T_cell = T_code * (1.0 + cvT * rng.standard_normal(n_tot))
        timesB, actB = gen_train(np.full(n_tot, N_FIXED), T_cell, n_tot, rng)
        # OLS 周期估计（归档口径）
        T_est = np.full(n_tot, np.nan)
        for i in range(n_tot):
            k = int(actB[i].sum())
            if k >= 3:
                T_est[i] = np.polyfit(np.arange(k), timesB[i, :k], 1)[0]
        ok = ~np.isnan(T_est)
        I_B_enc = mi_discrete(d[ok], qbin(T_est[ok], 12))
        m21B, mpuB = decoder(timesB, actB)
        I_B_dec = mi_discrete(d, qbin(m21B, 12) * 12 + qbin(mpuB, 12))
        B_rows.append(dict(cvT=cvT, I_enc=I_B_enc, I_dec=I_B_dec,
                           frac_estimable=float(ok.mean())))
        log(f"  臂B FM @CV(T)={cvT:.3f}: I(D;T_hat)={I_B_enc:.4f} bits, "
            f"过解码器 I={I_B_dec:.4f}")

    # ---- 臂C：幅度调制（AM） ----
    A_code = A0 * (1.0 + AMP_SLOPE * np.log(D))
    SK_LIST = [0.0, 0.02, 0.05, 0.075, 0.1, 0.2, 0.3, 0.45, 0.6]  # 加密低端网格以定位 AM<->计数交叉点
    C_rows = []
    for sk in SK_LIST:
        k = np.exp(rng.standard_normal(n_tot) * sk - 0.5 * sk * sk)
        # 过解码器：幅度串进积分器；逐脉冲内生噪声 AMP_INTR（与编码器读出同口径）
        timesC, actC = gen_train(np.full(n_tot, N_FIXED), np.full(n_tot, T_PERIOD), n_tot, rng)
        ampsC = np.where(actC, A_code[:, None] * k[:, None] *
                         (1.0 + AMP_INTR * rng.standard_normal(timesC.shape)), np.nan)
        # 编码器出口读出：同一组逐脉冲幅度的细胞均值（N_FIXED 脉冲取平均，口径一致）
        amp_meas = np.nanmean(ampsC, axis=1) + ETA_OBS * rng.standard_normal(n_tot)
        okc = ~np.isnan(amp_meas)  # 极少数首脉冲晚于观测窗的细胞（~0.02%），剔除以免 NaN 毒化分箱
        I_C_enc = mi_discrete(d[okc], qbin(amp_meas[okc], 16))
        with np.errstate(invalid="ignore"):
            kernel21 = Q21 * TAU21 * (1.0 - np.exp(-(W - timesC) / TAU21))
            m21C = np.nansum(ampsC * kernel21, axis=1)
            pidx = np.arange(N_MAX)[None, :]
            pm = actC & (pidx >= PUMA_NTHR - 1)
            kernelPU = QPU * TAUPU * (1.0 - np.exp(-(W - timesC) / TAUPU))
            mpuC = np.nansum(np.where(pm, ampsC, 0.0) * kernelPU, axis=1)
        I_C_dec = mi_discrete(d, qbin(m21C, 12) * 12 + qbin(mpuC, 12))
        assert I_C_dec <= H_D + 0.05, f"data-processing violation at sigma_k={sk}"
        C_rows.append(dict(sigma_k=sk, I_enc=I_C_enc, I_dec=I_C_dec))
        log(f"  臂C AM @sigma_k={sk:.2f}: I(D;A_hat)={I_C_enc:.4f} bits, "
            f"过解码器 I={I_C_dec:.4f}")

    # ---- 归档交叉验证：臂C 应复现代码81 假想幅度臂口径 ----
    # 代码81：sigma_k=0.6 时 I_amp=0.096 bits（假想臂，幅度份额 5.7%）
    return dict(armA=dict(I_enc=I_A_enc, I_dec=I_A_dec,
                          mean_pulses=float(n_pulse_A.mean())),
                armB_fm=B_rows, armC_am=C_rows,
                note_budget=f"臂B/C 固定 N={N_FIXED} 脉冲（对 FM/AM 有利的预算偏置）")

# ---------------- 主流程 ----------------

def main():
    log = lambda *a: print(*a, flush=True)
    log(f"[代码83] 种子 {SEED}，H(D)={H_D:.3f} bits")
    A = part_A()
    B = part_B()

    # ---- 汇总判词逻辑 ----
    c50 = [r for r in A["c_sweep"] if r["c"] == int(C_CLUSTER)][0]
    opt = A["optimal_quantizer"]
    # B 竞赛在"天然散布"点的对比：FM @ CV(T)=0.016（归档钉死值）；AM @ sigma_k=0.3（S5.6 细胞蛋白变异）
    fm_native = [r for r in B["armB_fm"] if abs(r["cvT"] - 0.016) < 1e-9][0]
    am_native = [r for r in B["armC_am"] if abs(r["sigma_k"] - 0.3) < 1e-9][0]
    am_06 = [r for r in B["armC_am"] if abs(r["sigma_k"] - 0.6) < 1e-9][0]
    verdict = dict(
        A_count_enc=B["armA"]["I_enc"], A_count_dec=B["armA"]["I_dec"],
        FM_native_enc=fm_native["I_enc"], FM_native_dec=fm_native["I_dec"],
        AM_native_enc=am_native["I_enc"], AM_native_dec=am_native["I_dec"],
        AM_sk06_enc=am_06["I_enc"],
        count_wins_enc=B["armA"]["I_enc"] > max(fm_native["I_enc"], am_native["I_enc"]),
        count_wins_dec=B["armA"]["I_dec"] > max(fm_native["I_dec"], am_native["I_dec"]),
    )
    log(f"  [判词] 天然散布点: 计数 {verdict['A_count_enc']:.4f} vs "
        f"FM {verdict['FM_native_enc']:.4f} vs AM {verdict['AM_native_enc']:.4f} bits "
        f"(编码器出口)；过解码器 {verdict['A_count_dec']:.4f} / "
        f"{verdict['FM_native_dec']:.4f} / {verdict['AM_native_dec']:.4f}")
    log(f"  [判词] 计数臂是否双赢: enc={verdict['count_wins_enc']} dec={verdict['count_wins_dec']}")

    double_record = [
        "Part A 全部用解析 Poisson 信道（非 MC），与代码81 的 MC 口径（I(D;N_DSB)=2.526）并列报告；两者差异为估计口径差异",
        "最优量化器为精确 DP（最小化 sum p(b) H(D|b)），M=13 级对齐 N_MAX=12+1；它给出的是'脉冲帽约束下的信息上限'，不代表生物可实现性",
        "修复稀释臂 phi 为自由参数（DSB 在触发前被修复的存活比例），文献半衰期范围未锁死，按扫描报告，结论只在'phi 越小损失越大'的单调层面引用",
        "臂B/C 固定 N=5 脉冲：臂A 全体平均脉冲约 2.5，该选择对 FM/AM 臂有利（更多脉冲=更好的周期/幅度估计），属保守（不利于计数臂）的公平性偏置",
        "FM 臂 T(D) 映射（gamma=0.15，2 倍动态范围）为任意选择；扫描显示结论对该斜率的方向不敏感的部分在 JSON 中，斜率本身未优化",
        "AM 臂直接复用代码81 假想幅度臂口径（斜率 0.15/ln、内生 CV 0.15、LogNormal 增益），sigma_k=0.6 行应与代码81 的 I_amp=0.096 bits 同量级",
        "解码级三臂共用同一对泄漏积分器（tau=10h/4h）：这本身偏袒计数臂，因为积分器天然是计数器；该偏置正是论点的一部分（已知解码器硬件是积分器），但登记为结构性偏置",
        "互信息为分箱 plug-in 估计（12x12 / 16 分位分箱），绝对值有偏、趋势可靠（代码50/81 同口径声明）",
    ]

    results = dict(
        meta=dict(script="代码83_为什么数字化_变分竞赛.py", seed=SEED, date="2026-09-25",
                  question="Q-A: L1 损失 0.872 bits 的物理不可约部分 vs 设计部分；"
                           "Q-B: 同一上游噪声下计数/FM/AM 三编码器结构的变分竞赛",
                  constants=dict(K0_DSB=K0_DSB, C_CLUSTER=C_CLUSTER, N_MAX=N_MAX,
                                 T_PERIOD=T_PERIOD, JIT_IPI=JIT_IPI, W=W,
                                 N_FIXED=N_FIXED, GAMMA_T=GAMMA_T,
                                 AMP_SLOPE=AMP_SLOPE, AMP_INTR=AMP_INTR)),
        H_D=H_D, part_A=A, part_B=B, verdict=verdict,
        double_record=double_record,
    )
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    def _js(o):
        if isinstance(o, np.bool_):
            return bool(o)
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            return float(o)
        return str(o)
    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2, default=_js)
    log("JSON 已写出:", OUT_JSON)

    # ---------------- 六面体图 ----------------
    fig = plt.figure(figsize=(17, 10.5))
    gs = fig.add_gridspec(2, 3, hspace=0.44, wspace=0.32)

    # (a) 聚类常数扫描 + 最优量化器上限
    axa = fig.add_subplot(gs[0, 0])
    cs = [r["c"] for r in A["c_sweep"]]
    Is = [r["I"] for r in A["c_sweep"]]
    axa.plot(cs, Is, "o-", color="#33527a", label="ceil($N_{DSB}$/c) clustering")
    axa.axhline(A["I_DSB_exact"], color="#0b6b3a", ls="--", lw=1.2,
                label=f"perfect sensing {A['I_DSB_exact']:.3f} bits")
    axa.axhline(opt["I"], color="#8c1d18", ls=":", lw=1.4,
                label=f"optimal 13-level quantizer {opt['I']:.3f}")
    axa.axvline(50, color="#c05640", ls="-", lw=0.9, alpha=0.6)
    axa.text(40, 0.62, "c=50\n(code 81)", fontsize=8, color="#c05640", ha="right")
    axa.set_xscale("log"); axa.set_xlabel("clustering constant c (DSB/trigger)")
    axa.set_ylabel("I(D; $N_{trig}$)  bits"); axa.legend(fontsize=8, loc="lower right")
    axa.set_title("(a) L1 mapping: is c=50 information-optimal?", fontsize=10.5)

    # (b) 修复稀释扫描
    axb = fig.add_subplot(gs[0, 1])
    phis = [r["phi"] for r in A["thinning_sweep"]]
    axb.plot(phis, [r["I_perfect_sense"] for r in A["thinning_sweep"]], "o-",
             color="#0b6b3a", label="perfect sensing after repair")
    axb.plot(phis, [r["I_optimal_quant"] for r in A["thinning_sweep"]], "s-",
             color="#8c1d18", label="+ optimal 13-level quantizer")
    axb.plot(phis, [r["I_c50"] for r in A["thinning_sweep"]], "^--",
             color="#c05640", label="+ c=50 clustering (code 81)")
    axb.axhline(A["I_DSB_exact"], color="#555", ls=":", lw=1.0)
    axb.set_xlabel("DSB survival fraction $\\varphi$ at sensing time")
    axb.set_ylabel("I(D; L1 output)  bits"); axb.legend(fontsize=8, loc="lower left")
    axb.set_title("(b) Repair thinning: the physical part of the L1 loss", fontsize=10.5)

    # (c) L1 损失分解条形
    axc = fig.add_subplot(gs[0, 2])
    phi_ref = 0.25  # 参考稀释点（双录：任意选择，只作分解示意）
    row = [r for r in A["thinning_sweep"] if abs(r["phi"] - phi_ref) < 1e-9][0]
    parts = [H_D - A["I_DSB_exact"], row["loss_thinning"],
             row["loss_total_optimal"] - row["loss_thinning"],
             row["loss_total_c50"] - row["loss_total_optimal"], row["I_c50"]]
    labels = ["Poisson\nfloor (L0)", f"repair\nthinning\n$\\varphi$={phi_ref}",
              "quantization\n(unavoidable\nunder cap)", "c=50\nsuboptimality",
              "surviving\nI(D;N_trig)"]
    colors = ["#33527a", "#d9863d", "#c05640", "#8c1d18", "#0b6b3a"]
    run = 0.0
    for i, (p_, l_, co_) in enumerate(zip(parts, labels, colors)):
        axc.bar(i, p_, bottom=run, color=co_, width=0.66)
        axc.text(i, run + p_ / 2, f"{p_:.2f}", ha="center", va="center", fontsize=8.5,
                 color="white" if i != 4 else "white", fontweight="bold")
        run += p_
    axc.axhline(H_D, color="#555", ls=":", lw=1.0)
    axc.text(0.0, H_D + 0.03, f"H(D)={H_D:.2f}", fontsize=8, ha="left")
    axc.set_ylim(0, H_D * 1.12)
    axc.set_xticks(range(5)); axc.set_xticklabels(labels, fontsize=7.2)
    axc.set_ylabel("bits")
    axc.set_title("(c) L1 loss decomposition (illustrative $\\varphi$=0.25)", fontsize=10.5)

    # (d) 编码器出口竞赛
    axd = fig.add_subplot(gs[1, 0])
    axd.plot([r["sigma_k"] for r in B["armC_am"]], [r["I_enc"] for r in B["armC_am"]],
             "s-", color="#c05640", label="AM arm vs gain $\\sigma_k$")
    axd.plot([r["cvT"] for r in B["armB_fm"]], [r["I_enc"] for r in B["armB_fm"]],
             "^-", color="#d9863d", label="FM arm vs CV(T)")
    axd.axhline(B["armA"]["I_enc"], color="#0b6b3a", ls="-", lw=1.6,
                label=f"counting arm (native) {B['armA']['I_enc']:.3f}")
    axd.axvline(0.016, color="#d9863d", ls=":", lw=0.9)
    axd.axvline(0.3, color="#c05640", ls=":", lw=0.9)
    axd.text(0.017, 0.05, "archived\nCV(T)=0.016", fontsize=7.5, color="#d9863d")
    axd.text(0.305, 0.05, "cellular $\\sigma_k$=0.3", fontsize=7.5, color="#c05640")
    axd.set_xlabel("native dispersion of the coding statistic")
    axd.set_ylabel("I(D; statistic)  bits"); axd.legend(fontsize=8, loc="upper right")
    axd.set_title("(d) Encoder-output competition under native noise", fontsize=10.5)

    # (e) 过解码器竞赛
    axe = fig.add_subplot(gs[1, 1])
    axe.plot([r["sigma_k"] for r in B["armC_am"]], [r["I_dec"] for r in B["armC_am"]],
             "s-", color="#c05640", label="AM arm")
    axe.plot([r["cvT"] for r in B["armB_fm"]], [r["I_dec"] for r in B["armB_fm"]],
             "^-", color="#d9863d", label="FM arm")
    axe.axhline(B["armA"]["I_dec"], color="#0b6b3a", ls="-", lw=1.6,
                label=f"counting arm {B['armA']['I_dec']:.3f}")
    axe.set_xlabel("native dispersion of the coding statistic")
    axe.set_ylabel("I(D; mRNA joint)  bits"); axe.legend(fontsize=8)
    axe.set_title("(e) Same competition through the leaky-integrator decoder", fontsize=10.5)

    # (f) 摘要
    axs = fig.add_subplot(gs[1, 2]); axs.axis("off")
    v = verdict
    summary = (
        f"Code 83 verdict summary (seed {SEED})\n\n"
        f"Q-A  L1 loss decomposition:\n"
        f"  I(D;N_DSB) exact = {A['I_DSB_exact']:.3f} bits\n"
        f"  Poisson floor (L0) = {A['L0_poisson_floor']:.3f} bits\n"
        f"  c=50 clustering I = {c50['I']:.3f} (loss {c50['loss']:.3f})\n"
        f"  optimal 13-level I = {opt['I']:.3f} (loss {opt['loss']:.3f})\n"
        f"  -> c=50 suboptimality = {row['loss_total_c50'] - row['loss_total_optimal']:.3f} bits\n\n"
        f"Q-B  encoder race @ native noise:\n"
        f"  counting: enc {v['A_count_enc']:.3f} / dec {v['A_count_dec']:.3f} bits\n"
        f"  FM @CV(T)=0.016: enc {v['FM_native_enc']:.3f} / dec {v['FM_native_dec']:.3f}\n"
        f"  AM @$\\sigma_k$=0.3: enc {v['AM_native_enc']:.3f} / dec {v['AM_native_dec']:.3f}\n"
        f"  AM @$\\sigma_k$=0.6: enc {v['AM_sk06_enc']:.3f} (code81: 0.096)\n\n"
        f"  counting wins enc: {v['count_wins_enc']}\n"
        f"  counting wins dec: {v['count_wins_dec']}\n"
        f"  budget note: FM/AM arms use N=5 fixed pulses\n"
        f"  (bias against the counting arm, by design)"
    )
    axs.text(0.02, 0.98, summary, va="top", fontsize=8.4, family="monospace",
             transform=axs.transAxes)
    fig.suptitle("Code 83: why digital, why noise upstream — L1 loss decomposition + encoder variational race",
                 fontsize=12.5)
    fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
    fig.savefig(OUT_SVG, bbox_inches="tight")
    log("图已写出:", OUT_PNG, "和", OUT_SVG)
    log("[代码83] 完成。")


if __name__ == "__main__":
    main()
