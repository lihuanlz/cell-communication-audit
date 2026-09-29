# -*- coding: utf-8 -*-
"""
代码81：p53 全通链路级联模型与逐接口审计（L0 剂量 -> L4 命运）
=================================================================
日期：2026-09-23 ｜ 种子固定：20260923 ｜ 全程确定性

谱系：
  - 设计文档：03_细胞线3/文稿/设计_p53全通链路级联模型_v01.md
  - 编码器层已封卷：判词卡_P53-Fisher审计_2026-09-18.md / 封卷卡_p53信号线_2026-09-18.md
    （代码51r3 / 代码51r4A 归档参数）
  - 玩具级联自检：代码50（02_细胞线2，v0.2 双臂版）
  - 框架：Nature_SI_v05 S1.3（幸存集）/ S1.5（三腿判据）/ S1.7（信息损失恒等式）/ S5.6（数字化代价）

归档锚（逐字引用，重叠处必须一致）：
  - 时间类 Fisher 份额 94.2%，幅度 5.5%
  - I(D;N) = 3.07 / 3.17 bits（96.8%）
  - 增益腐蚀：sigma_k 0->0.6，CV(A) 0.002->0.639，CV(T) 钉死 0.016
  - sigma* = 0.25 -> 细胞内 IPI CV = 0.316 ≈ 文献 0.30；基线 sigma=0.06 -> IPI CV 0.085
  - mu_d ≈ 1.0 T -> Var(t1)/Var(T) = 5.7 ≈ 文献 5.8

模型层定义：
  L0 损伤层：N_DSB ~ Poisson(35 * D)          （计数产出，类别 I）
  L1 传感层：N_trig = ceil(N_DSB / 50)，cap 12 （registered gap，最简聚类映射，假设）
  L2 编码层：T = 5.5 h 固定常数；A 固定（CV 0.002）；N = N_trig；
             首脉冲延迟 t1 ~ Exp(1.0 T)；区间抖动 CV 0.085（归档 sigma=0.06 行）
  L3 解码层：mRNA 泄漏积分；p21 型积分器 tau=10 h；PUMA 型累积器 tau=4 h，第 3 脉冲起才产出
  L4 命运层：阈值于 PUMA/p21 积分比与绝对水平 -> {survival, arrest, apoptosis}
"""

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
try:
    from daimon_runtime import setup_plot
    setup_plot()
except Exception:
    pass

SEED = 20260923
ROOT = Path(__file__).resolve().parent.parent
OUT_JSON = ROOT / "结果" / "代码81_级联审计_结果.json"
OUT_PNG = ROOT / "结果" / "代码81_级联审计_五面体.png"

# ---------------- 固定常数（全部显式，不拟合） ----------------
DOSES = np.array([0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0])   # Gy，等权 -> H(D)=log2(7)
N_CELLS = 4000                                            # 每剂量细胞数
K0_DSB = 35.0         # DSB/Gy（文献锚）
C_CLUSTER = 50.0      # DSB/触发（L1 聚类常数，假设：10 Gy -> 7 脉冲，对齐 0-7 脉冲/24h 文献范围）
N_MAX = 12            # 触发上限
T_PERIOD = 5.5        # h，候选结构常数（类别 II），归档固定值
MU_D = 1.0 * T_PERIOD # 首脉冲延迟均值（归档 A3 标定：mu_d≈T 时比值 5.7≈文献 5.8）
JIT_IPI = 0.085       # 逐区间抖动 CV（归档 r4A sigma=0.06 行 IPI CV = 0.085）
A0 = 1.0              # 固定幅度（不载剂量信息）
A_CV = 0.002          # 幅度内生 CV（归档 CV(A)=0.002 @ sigma_k=0）
W = 48.0              # 观测窗 h
Q21, TAU21 = 100.0, 10.0        # p21 型积分器：每脉冲量子产出、泄漏时间常数
QPU, TAUPU, PUMA_NTHR = 100.0, 4.0, 3   # PUMA 型累积器：第 3 个脉冲起才产出（高阈值）
TH_PUMA, R0_FATE, TH_P21 = 800.0, 0.15, 1500.0  # L4 命运阈值（固定常数）
B_GAIN = 0.3          # 仿射观测基底 b（代码50 同口径）
ETA_OBS = 0.02        # 读出加性噪声（相对 A0）
N_BOOT = 200          # bootstrap 次数
ALPHA_E = 7.0         # S5.6 事件可靠性常数

# ---------------- 互信息工具 ----------------

def mi_discrete(x, y):
    """plug-in 离散互信息（bits），x,y 为整数标签数组。"""
    x = np.asarray(x); y = np.asarray(y)
    Kx = int(x.max()) + 1; Ky = int(y.max()) + 1
    N = len(x)
    cxy = np.bincount(x * Ky + y, minlength=Kx * Ky).reshape(Kx, Ky) / N
    cx = cxy.sum(1, keepdims=True); cy = cxy.sum(0, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        term = cxy * np.log2(cxy / (cx @ cy))
    return float(np.nansum(term))


def mi_mm(x, y):
    """Miller-Madow 偏差校正后的互信息（bits），下限 0。"""
    x = np.asarray(x); y = np.asarray(y)
    rx = len(np.unique(x)); ry = len(np.unique(y)); N = len(x)
    bias = (rx - 1) * (ry - 1) / (2 * N * np.log(2))
    return max(0.0, mi_discrete(x, y) - bias)


def qbin(v, nbin):
    """分位数分箱 -> 整数标签。"""
    v = np.asarray(v, dtype=float)
    edges = np.unique(np.quantile(v, np.linspace(0, 1, nbin + 1)))
    if len(edges) <= 2:
        return np.zeros(len(v), dtype=int)
    return np.clip(np.digitize(v, edges[1:-1]), 0, len(edges) - 2)


def boot_ci(x, y, n_boot=N_BOOT, rng=None):
    """bootstrap 百分位 CI（plug-in MI）。"""
    rng = rng or np.random.default_rng(0)
    N = len(x); vals = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, N, N)
        vals[b] = mi_discrete(x[idx], y[idx])
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))

# ---------------- 级联仿真 ----------------

def simulate(dose, n, rng, cvT_extra=0.0):
    """单剂量群体仿真。返回每层输出字典。

    cvT_extra：逐细胞周期散布（固定常数检验用，审计 4）。
    """
    # L0：损伤计数
    n_dsb = rng.poisson(K0_DSB * dose, n)
    # L1：聚类触发映射（registered gap，最简形式，假设）
    n_trig = np.minimum(np.ceil(n_dsb / C_CLUSTER).astype(int), N_MAX)
    # L2：脉冲串（事件时间）
    t1 = rng.exponential(MU_D, n)
    T_cell = T_PERIOD * (1.0 + cvT_extra * rng.standard_normal(n))  # 逐细胞周期散布
    jit = JIT_IPI
    times = np.full((n, N_MAX), np.nan)
    tcur = t1.copy()
    for j in range(N_MAX):
        times[:, j] = tcur
        interval = T_cell * (1.0 + jit * rng.standard_normal(n))
        interval = np.maximum(interval, 0.3 * T_PERIOD)
        tcur = tcur + interval
    in_win = times <= W
    pulse_idx = np.arange(N_MAX)[None, :]
    active = (pulse_idx < n_trig[:, None]) & in_win
    n_pulse = active.sum(1)
    # 幅度：固定 + 内生 CV 0.002
    amps = np.where(active, A0 * (1.0 + A_CV * rng.standard_normal(times.shape)), np.nan)
    # L3：泄漏积分解码
    with np.errstate(invalid="ignore"):
        integ = lambda q, tau: np.nansum(
            np.where(active, q * tau * (1.0 - np.exp(-(W - times) / tau)), 0.0), axis=1)
    m_p21 = integ(Q21, TAU21)
    puma_mask = active & (pulse_idx >= PUMA_NTHR - 1)
    m_puma = np.nansum(np.where(puma_mask,
                                QPU * TAUPU * (1.0 - np.exp(-(W - times) / TAUPU)), 0.0), axis=1)
    # L4：命运 = 计数比与水平的阈值
    ratio = m_puma / (m_p21 + 1.0)
    fate = np.where((m_puma >= TH_PUMA) & (ratio >= R0_FATE), 2,
                    np.where(m_p21 >= TH_P21, 1, 0))  # 0 survival, 1 arrest, 2 apoptosis
    # 逐细胞周期估计（事件时间对序号 OLS 斜率 -> CV(T) 口径与归档一致）
    T_est = np.full(n, np.nan)
    for i in range(n):
        k = int(n_pulse[i])
        if k >= 3:
            xs = np.arange(k)
            T_est[i] = np.polyfit(xs, times[i, :k], 1)[0]
    return dict(n_dsb=n_dsb, n_trig=n_trig, n_pulse=n_pulse, t1=t1,
                times=times, amps=amps, m_p21=m_p21, m_puma=m_puma,
                fate=fate, T_est=T_est)


def run_ensemble(cvT_extra=0.0, seed=SEED):
    rng = np.random.default_rng(seed)
    parts = []
    for di, D in enumerate(DOSES):
        p = simulate(D, N_CELLS, rng, cvT_extra=cvT_extra)
        p["dose_idx"] = np.full(N_CELLS, di)
        parts.append(p)
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}

# ---------------- 主流程 ----------------

def main():
    log = lambda *a: print(*a, flush=True)
    log("[代码81] 全链路级联仿真启动，种子", SEED)
    E = run_ensemble()
    d = E["dose_idx"]
    H_D = np.log2(len(DOSES))
    log(f"H(D) = {H_D:.3f} bits, 细胞总数 {len(d)}")

    # ---- 命运分布检查 ----
    fate_names = ["survival", "arrest", "apoptosis"]
    fate_tab = {}
    for di, D in enumerate(DOSES):
        m = d == di
        fr = [float(np.mean(E["fate"][m] == k)) for k in range(3)]
        fate_tab[f"{D}Gy"] = dict(zip(fate_names, [round(f, 4) for f in fr]))
        log(f"  D={D:>5} Gy  fate 分布 S/A/P = {[round(f,3) for f in fr]}, "
            f"<N_DSB>={E['n_dsb'][m].mean():.1f}, <N_pulse>={E['n_pulse'][m].mean():.2f}")

    # ---- 审计 2：信息损失恒等式 H(D) = I(D;Y5) + sum L_i ----
    y1 = E["n_dsb"]                       # L0 输出
    y2 = E["n_trig"]                      # L1 输出
    y3 = E["n_pulse"] * 8 + qbin(E["t1"], 8)          # L2 输出（计数 + 首脉冲延迟分箱）
    y4 = qbin(E["m_p21"], 12) * 12 + qbin(E["m_puma"], 12)  # L3 输出（联合分箱）
    y5 = E["fate"]                        # L4 输出
    labels = ["Y1=N_DSB (L0)", "Y2=N_trig (L1)", "Y3=pulse train (L2)",
              "Y4=mRNA joint (L3)", "Y5=fate (L4)"]
    ys = [y1, y2, y3, y4, y5]
    rng_b = np.random.default_rng(SEED + 1)
    I_chain, I_chain_mm, I_ci = [], [], []
    for lab, y in zip(labels, ys):
        I0 = mi_discrete(d, y)
        lo, hi = boot_ci(d, y, rng=rng_b)
        I_chain.append(I0); I_chain_mm.append(mi_mm(d, y)); I_ci.append([lo, hi])
        log(f"  I(D;{lab}) = {I0:.4f} bits  MM={mi_mm(d,y):.4f}  CI95=[{lo:.4f},{hi:.4f}]")
    I_seq = [H_D] + I_chain               # I(D;Y0)=H(D)
    losses = [I_seq[i] - I_seq[i + 1] for i in range(5)]
    closure = I_chain[-1] + sum(losses)
    weakest = int(np.argmax(losses))
    weakest_name = ["L0 damage", "L1 sensing(gap)", "L2 encoder", "L3 decoder", "L4 fate"][weakest]
    log(f"  恒等式闭合：I(D;Y5)+sum(L_i) = {closure:.6f} vs H(D) = {H_D:.6f}, "
        f"残差 {closure - H_D:.2e}")
    log(f"  逐层损失 L0..L4 = {[round(float(l),4) for l in losses]}, 最弱层 = {weakest_name}")

    # ---- 审计 3：增益散布测试（计数链 vs 幅度链） ----
    SIG_K = [0.0, 0.1, 0.2, 0.3, 0.45, 0.6]
    rng_g = np.random.default_rng(SEED + 2)
    gain_sweep = []
    amps = E["amps"]; times = E["times"]; n_pulse = E["n_pulse"]
    # 假想幅度臂：A(D) = A0*(1+0.15*ln(D/1Gy))，幅度跨剂量范围约翻倍（代码50 臂B 口径）
    amp_code = A0 * (1.0 + 0.15 * np.log(DOSES[d]))[:, None] * np.ones_like(amps)
    amp_intr = 0.15  # 假想幅度臂内生 CV（代码50 臂A 口径 1.0±0.15）
    for sk in SIG_K:
        k = np.exp(rng_g.standard_normal(len(d)) * sk - 0.5 * sk * sk)  # LogNormal(0,sk)，中位 1
        y_obs = k[:, None] * amps + B_GAIN + ETA_OBS * rng_g.standard_normal(amps.shape)
        # 计数检测：逐细胞仿射不变阈值（无脉冲期基底估计 + 0.5*动态范围）
        base = B_GAIN + ETA_OBS * rng_g.standard_normal((len(d), 8))
        b_est = np.median(base, axis=1)
        has = n_pulse > 0
        y_max = np.full(len(d), np.nan)
        amp_meas = np.full(len(d), np.nan)
        y_max[has] = np.nanmax(y_obs[has], axis=1)
        amp_meas[has] = k[has] * np.nanmean(amps[has], axis=1)
        thr = b_est + 0.5 * (y_max - b_est)
        n_det = np.where(has, np.nansum(y_obs > thr[:, None], axis=1), 0).astype(int)
        # 幅度链（假想臂）：测幅度 = k * A_code * (1+0.15 z)
        amp_hyp = k * amp_code[:, 0] * (1.0 + amp_intr * rng_g.standard_normal(len(d)))
        amp_hyp = amp_hyp + ETA_OBS * rng_g.standard_normal(len(d))
        I_cnt = mi_discrete(d, n_det)
        I_amp = mi_discrete(d, qbin(amp_hyp, 16))
        # 编码器口径复现：CV(A_meas) 与 CV(T_est)
        real = ~np.isnan(amp_meas)
        cv_A = float(np.std(amp_meas[real]) / np.mean(amp_meas[real]))
        tval = E["T_est"]; realt = ~np.isnan(tval)
        cv_T = float(np.std(tval[realt]) / np.mean(tval[realt]))
        # 归档标定子集：N_pulse >= 7（对齐归档 N≈9 的周期估计口径）
        calib = realt & (n_pulse >= 7)
        cv_T_cal = float(np.std(tval[calib]) / np.mean(tval[calib])) if calib.sum() > 50 else np.nan
        gain_sweep.append(dict(sigma_k=sk, I_count=I_cnt, I_amp_hyp=I_amp,
                               cv_A=cv_A, cv_T=cv_T, cv_T_calib=cv_T_cal,
                               frac_count_intact=float(np.mean(n_det == n_pulse))))
        log(f"  sigma_k={sk:.2f}: I(count)={I_cnt:.4f}  I(amp_hyp)={I_amp:.4f}  "
            f"CV(A)={cv_A:.3f}  CV(T)={cv_T:.4f}/{cv_T_cal:.4f}(N>=7)  "
            f"计数检出率={np.mean(n_det==n_pulse):.4f}")

    # ---- 审计 4：固定常数检验（CV(T) 散布扫描） ----
    CVT = [0.0, 0.016, 0.05, 0.1, 0.2]
    cvt_sweep = []
    for cv in CVT:
        Ec = run_ensemble(cvT_extra=cv, seed=SEED + 3)
        I_N = mi_discrete(Ec["dose_idx"], Ec["n_pulse"])
        I_F = mi_discrete(Ec["dose_idx"], Ec["fate"])
        cvt_sweep.append(dict(cvT=cv, I_count=I_N, I_fate=I_F))
        log(f"  CV(T)={cv:.3f}: I(D;N)={I_N:.4f}  I(D;fate)={I_F:.4f}")
    # 假想幅度码 vs 计数码的端到端份额（标定增益散布 sigma_k=0.3，S5.6 细胞蛋白变异 ~30%）
    rng_a = np.random.default_rng(SEED + 4)
    amp_share = {}
    for sk in [0.0, 0.3, 0.6]:
        k = np.exp(rng_a.standard_normal(len(d)) * sk - 0.5 * sk * sk)
        amp_hyp = k * amp_code[:, 0] * (1.0 + amp_intr * rng_a.standard_normal(len(d)))
        I_amp = mi_mm(d, qbin(amp_hyp, 16))
        I_cnt = mi_mm(d, E["n_pulse"])
        amp_share[f"sigma_k={sk}"] = dict(I_amp=I_amp, I_count=I_cnt,
                                          share=I_amp / (I_amp + I_cnt))
        log(f"  幅度份额 @sigma_k={sk}: I_amp={I_amp:.4f} vs I_cnt={I_cnt:.4f} "
            f"-> share={I_amp/(I_amp+I_cnt):.4f}")

    # ---- 审计 5：能耗核算（S5.6） ----
    bits_end = mi_mm(d, y5)
    D_levels = float(2 ** bits_end)     # 端到端可分辨级数（由传输比特反演）
    L_axis = np.arange(1, 9)
    E_A = 2.0 * L_axis**2 * D_levels**2    # 2 kT L^2 D^2
    E_B = ALPHA_E * L_axis * D_levels**2   # alpha kT L D^2
    L_star = ALPHA_E / 2.0
    L_cascade = 4
    energy = dict(bits_per_cell=bits_end, D_levels=D_levels,
                  E_A_at_L4=float(2 * L_cascade**2 * D_levels**2),
                  E_B_at_L4=float(ALPHA_E * L_cascade * D_levels**2),
                  ratio_EA_over_EB_at_L4=float(2 * L_cascade / ALPHA_E),
                  L_star=L_star,
                  E_A_series=E_A.tolist(), E_B_series=E_B.tolist(), L_axis=L_axis.tolist())
    log(f"  能耗：端到端 {bits_end:.4f} bits/细胞 -> D={D_levels:.2f} 级；"
        f"L=4 处 E_A/E_B = {2*L_cascade/ALPHA_E:.3f}，交叉 L*={L_star}")

    # ---- 审计 1：三腿判据 x 4 接口 ----
    g0, g6 = gain_sweep[0], gain_sweep[-1]
    legs = []
    def leg(name, legs_status, evidence, flagged=False):
        legs.append(dict(interface=name, legs=legs_status, evidence=evidence,
                         assumption_flagged=flagged))
        log(f"  接口 {name}: {legs_status} | {evidence}")
    leg("I0->1 剂量->DSB计数",
        dict(leg1_algebraic=True, leg2_symbolic=True, leg3_temporal=True),
        f"计数为类别I免校准统计量（S1.4）；I(D;N_DSB)={I_chain[0]:.3f} bits > 0；"
        f"DSB 形成时标（分钟）远小于剂量交付与下游时标（小时）")
    leg("I1->2 DSB->脉冲触发（registered gap）",
        dict(leg1_algebraic=True, leg2_symbolic=True, leg3_temporal=True),
        f"聚类映射计数->计数，代数幸存平凡成立（假设）；I(D;N_trig)={I_chain[1]:.3f} bits；"
        f"ATM 激活分钟级 << T=5.5h；本层无公开单细胞数据，三腿均标注为假设",
        flagged=True)
    leg("I2->3 脉冲串->mRNA（计数/事件时间）",
        dict(leg1_algebraic=True, leg2_symbolic=True, leg3_temporal=True),
        f"计数对仿射简并免疫：sigma_k 0->0.6 检出率 {g0['frac_count_intact']:.3f}->"
        f"{g6['frac_count_intact']:.3f}，I(count) {g0['I_count']:.3f}->{g6['I_count']:.3f} bits；"
        f"解码器 tau_p21=10h、tau_PUMA=4h 夹住脉冲间隔 5.5h 与窗 48h（tau>W 失忆，tau<<T 无积分）")
    leg("I2->3 脉冲串->mRNA（幅度，反例腿）",
        dict(leg1_algebraic=False, leg2_symbolic=False, leg3_temporal=None),
        f"幅度腿1 死亡：sigma_k 0->0.6，CV(A_meas) {g0['cv_A']:.3f}->{g6['cv_A']:.3f}；"
        f"假想幅度码份额 @sigma_k=0.6 = {amp_share['sigma_k=0.6']['share']:.3f}"
        f"（归档编码器层 Fisher 份额 5.5% 的端到端对应）")
    leg("I3->4 mRNA->命运",
        dict(leg1_algebraic=True, leg2_symbolic=True, leg3_temporal=True),
        f"命运=计数比阈值，比值消共模增益（S1.4 差分类幸存）；I(D;fate)={I_chain[4]:.3f} bits > 0；"
        f"决策时标 24-48h 与 mRNA 积分时标匹配")

    # ---- 归档一致性核对（重叠数字逐字对照） ----
    archived_check = dict(
        cv_A_at_sk0=dict(model=g0["cv_A"], archived=0.002,
                         match=abs(g0["cv_A"] - 0.002) < 0.005),
        cv_A_at_sk06=dict(model=g6["cv_A"], archived=0.639, note="定性复现（生成模型 LogNormal 增益）",
                          match=0.55 <= g6["cv_A"] <= 0.75),
        cv_T_pinned=dict(model_range=[min(g["cv_T"] for g in gain_sweep),
                                      max(g["cv_T"] for g in gain_sweep)],
                         model_calib_N7=[round(g["cv_T_calib"], 4) for g in gain_sweep],
                         archived=0.016,
                         note="OLS 周期估计口径下钉死在小值且对增益散布不变；N>=7 标定子集对齐归档 N≈9 口径",
                         match=max(g["cv_T"] for g in gain_sweep) < 0.06),
        count_immune=dict(I_count_range=[min(g["I_count"] for g in gain_sweep),
                                         max(g["I_count"] for g in gain_sweep)],
                          note="计数链信息对 sigma_k 平坦"),
        archived_anchor=dict(fisher_time_share=0.942, fisher_amp_share=0.055,
                             I_DN="3.07/3.17 bits (96.8%)", cv_T_pinned=0.016,
                             sigma_star=0.25, ipi_cv_star=0.316, var_ratio_lit=5.7))
    for key in ["cv_A_at_sk0", "cv_A_at_sk06", "cv_T_pinned"]:
        mv = archived_check[key].get("model", archived_check[key].get("model_range"))
        log(f"  归档核对 {key}: model={mv} match={archived_check[key]['match']}")

    # ---- 失败与假设双录 ----
    double_record = [
        "L1 传感层为 registered gap：无公开单细胞 ATM 活性数据，聚类常数 c=50 DSB/触发为假设（选取依据：10 Gy -> 7 脉冲，对齐文献 0-7 脉冲/24h 剂量范围）；恒等式中 L2 损失含该层不可知成分",
        "L3/L4 解码核（量子产出、tau、PUMA 第 3 脉冲阈值、命运阈值 TH_PUMA=800/R0=0.15/TH_P21=1500）为文献启发假设，非拟合",
        "互信息为分箱 plug-in 估计（mRNA 12x12 分位分箱），绝对值有偏、趋势可靠（代码50 同口径声明）；同时报 Miller-Madow 校正值",
        "p21 积分器 tau=10h > 脉冲间隔 5.5h：逐脉冲分辨在该接口死亡，幸存的是积分计数（与代码50 腿3 tau=8h 崩溃一致），积分统计量的时标为整串时长而非单脉冲间隔",
        "幅度臂为假想反事实臂（A(D) 对数斜率 0.15、内生 CV 0.15），用于给幅度链一个'有信息可死'的对照；真实 p53 幅度剂量不变（归档 |alpha|<=0.05）",
        "CV(T) 的归档值 0.016 取自 FHN 代表元 sigma=0.06 行；本生成模型用 OLS 周期估计复现同量级钉死值，数值一致性依赖估计口径，已显式声明",
    ]

    results = dict(
        meta=dict(script="代码81_p53全通链路级联审计.py", seed=SEED, date="2026-09-23",
                  doses=DOSES.tolist(), n_cells_per_dose=N_CELLS,
                  constants=dict(K0_DSB=K0_DSB, C_CLUSTER=C_CLUSTER, N_MAX=N_MAX,
                                 T_PERIOD=T_PERIOD, MU_D_over_T=MU_D / T_PERIOD,
                                 JIT_IPI=JIT_IPI, A0=A0, A_CV=A_CV, W=W,
                                 Q21=Q21, TAU21=TAU21, QPU=QPU, TAUPU=TAUPU,
                                 PUMA_NTHR=PUMA_NTHR, TH_PUMA=TH_PUMA,
                                 R0_FATE=R0_FATE, TH_P21=TH_P21, B_GAIN=B_GAIN,
                                 ETA_OBS=ETA_OBS, ALPHA_E=ALPHA_E)),
        fate_table=fate_tab,
        info_identity=dict(H_D=H_D, labels=labels, I_chain=I_chain,
                           I_chain_mm=I_chain_mm, I_ci95=I_ci,
                           losses=losses, closure_sum=closure,
                           closure_residual=closure - H_D,
                           weakest_layer=weakest_name,
                           weakest_loss=losses[weakest]),
        gain_sweep=gain_sweep,
        cvt_sweep=cvt_sweep,
        amp_share=amp_share,
        energy=energy,
        three_legs=legs,
        archived_check=archived_check,
        double_record=double_record,
    )
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)
    log("JSON 已写出:", OUT_JSON)

    # ---------------- 五面体图 ----------------
    fig = plt.figure(figsize=(17, 10.5))
    gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.30)

    # (a) 级联图 + 幸存统计量标签
    axa = fig.add_subplot(gs[0, 0])
    axa.axis("off")
    boxes = [("L0 damage\n$N_{DSB}$ ~ Pois(35D)", "count"),
             ("L1 ATM sensing\n$N_{trig}$=ceil($N_{DSB}$/50)\n[registered gap]", "count (assumed)"),
             ("L2 p53 encoder\nT=5.5h fixed, A fixed", "count + event times\n(amplitude dead)"),
             ("L3 promoter\nleaky integral", "count integral"),
             ("L4 fate\nratio threshold", "ordinal fate")]
    ys_pos = np.linspace(0.88, 0.12, 5)
    for (txt, surv), yp in zip(boxes, ys_pos):
        axa.add_patch(plt.Rectangle((0.06, yp - 0.075), 0.52, 0.15,
                                    fc="#eef3fb", ec="#33527a", lw=1.4,
                                    transform=axa.transAxes))
        axa.text(0.32, yp, txt, ha="center", va="center", fontsize=8.2,
                 transform=axa.transAxes)
        axa.text(0.64, yp, "survives:\n" + surv, ha="left", va="center",
                 fontsize=7.6, color="#0b6b3a", transform=axa.transAxes)
    for i in range(4):
        axa.annotate("", xy=(0.32, ys_pos[i + 1] + 0.078), xytext=(0.32, ys_pos[i] - 0.078),
                     arrowprops=dict(arrowstyle="-|>", color="#33527a", lw=1.6),
                     xycoords=axa.transAxes)
    axa.set_title("(a) p53 cascade: surviving statistic per interface", fontsize=10.5)

    # (b) 信息瀑布
    axb = fig.add_subplot(gs[0, 1])
    names = ["H(D)", "L0\ndamage", "L1\nsense(gap)", "L2\nencode", "L3\ndecode", "L4\nfate", "I(D;fate)"]
    starts = [0.0]
    heights = [H_D]
    run = H_D
    for l in losses:
        starts.append(run - l); heights.append(l); run -= l
    starts.append(0.0); heights.append(I_chain[-1])
    colors = ["#33527a"] + ["#c05640"] * 4 + ["#d9863d"] + ["#0b6b3a"]
    colors[1 + weakest] = "#8c1d18"
    axb.bar(range(7), heights, bottom=starts, color=colors, width=0.68)
    for i, (s, h) in enumerate(zip(starts, heights)):
        axb.text(i, s + h + 0.04, f"{h:.2f}", ha="center", fontsize=8.5)
    axb.annotate(f"weakest: {weakest_name}\nL={losses[weakest]:.2f} bits",
                 xy=(1 + weakest, starts[1 + weakest] + 0.45 * heights[1 + weakest]),
                 xytext=(1 + weakest + 0.7, H_D * 1.16), fontsize=8.5, color="#8c1d18",
                 arrowprops=dict(arrowstyle="->", color="#8c1d18"))
    axb.set_xticks(range(7)); axb.set_xticklabels(names, fontsize=7.2)
    axb.set_ylabel("bits"); axb.set_ylim(0, H_D * 1.32)
    axb.set_title("(b) Information waterfall: $H(D)=I(D;Y_4)+\\sum L_i$", fontsize=10.5)

    # (c) 增益散布不对称
    axc = fig.add_subplot(gs[0, 2])
    sks = [g["sigma_k"] for g in gain_sweep]
    axc.plot(sks, [g["I_count"] for g in gain_sweep], "o-", color="#0b6b3a",
             label="counting chain I(D;N)")
    axc.plot(sks, [g["I_amp_hyp"] for g in gain_sweep], "s-", color="#c05640",
             label="hypothetical amplitude chain")
    axc.set_xlabel("gain dispersion $\\sigma_k$"); axc.set_ylabel("bits")
    axc2 = axc.twinx()
    axc2.plot(sks, [g["cv_A"] for g in gain_sweep], "^--", color="#8c1d18", label="CV(A)")
    axc2.plot(sks, [g["cv_T"] for g in gain_sweep], "v--", color="#33527a", label="CV(T) pop.")
    axc2.plot(sks, [g["cv_T_calib"] for g in gain_sweep], "v-", color="#5a7fb5",
              label="CV(T) N>=7 calib.")
    axc2.set_ylabel("CV", color="#555")
    axc2.axhline(0.639, color="#8c1d18", ls=":", lw=0.9)
    axc2.axhline(0.016, color="#33527a", ls=":", lw=0.9)
    axc2.set_ylim(-0.09, 0.75)
    axc2.text(0.22, 0.605, "archived 0.639", fontsize=7.5, color="#8c1d18")
    axc2.text(0.28, -0.055, "archived CV(T)=0.016", fontsize=7.5, color="#33527a")
    h1, l1 = axc.get_legend_handles_labels(); h2, l2 = axc2.get_legend_handles_labels()
    axc.legend(h1 + h2, l1 + l2, fontsize=7.5, loc="center left")
    axc.set_title("(c) Gain-dispersion asymmetry: count survives, amplitude dies", fontsize=10.5)

    # (d) 固定常数鲁棒性
    axd = fig.add_subplot(gs[1, 0])
    cvs = [c["cvT"] for c in cvt_sweep]
    axd.plot(cvs, [c["I_count"] for c in cvt_sweep], "o-", color="#0b6b3a",
             label="I(D; pulse count)")
    axd.plot(cvs, [c["I_fate"] for c in cvt_sweep], "s-", color="#33527a",
             label="I(D; fate) end-to-end")
    axd.axvline(0.016, color="#8c1d18", ls=":", lw=1.0)
    axd.text(0.017, axd.get_ylim()[0] + 0.05, "archived\nCV(T)=0.016", fontsize=7.5,
             color="#8c1d18")
    sh = amp_share["sigma_k=0.6"]
    axd.text(0.55, 0.18, f"hypothetical amplitude share @ $\\sigma_k$=0.6:\n"
                         f"{sh['share']*100:.1f}% (archived encoder share 5.5%)",
             transform=axd.transAxes, fontsize=8, color="#c05640",
             bbox=dict(fc="#fdf1ec", ec="#c05640", lw=0.8))
    axd.set_xlabel("per-cell CV(T) dispersion"); axd.set_ylabel("bits")
    axd.legend(fontsize=8, loc="upper right")
    axd.set_title("(d) Fixed-constant robustness: information vs CV(T)", fontsize=10.5)

    # (e) 能耗对比
    axe = fig.add_subplot(gs[1, 1])
    axe.plot(L_axis, E_A, "o-", color="#c05640", label="analogue $E_A=2kTL^2D^2$")
    axe.plot(L_axis, E_B, "s-", color="#0b6b3a", label="counting $E_B=\\alpha kTLD^2$, $\\alpha$=7")
    axe.axvline(L_star, color="#555", ls="--", lw=1.0)
    axe.text(L_star + 0.08, E_A[-1] * 0.75, f"crossover $L^*=\\alpha/2$={L_star}", fontsize=8.5)
    axe.axvline(L_cascade, color="#33527a", ls=":", lw=1.2)
    axe.text(L_cascade + 0.08, E_A[0] * 1.6, "p53 cascade L=4\n(counting cheaper by "
             f"{(2*L_cascade/ALPHA_E - 1)*100:.0f}%)", fontsize=8.5, color="#33527a")
    axe.set_xlabel("cascade depth L (transmission stages)")
    axe.set_ylabel(f"energy (kT units, D={D_levels:.2f} levels)")
    axe.legend(fontsize=8.5, loc="upper left")
    axe.set_title("(e) Energy: counting vs analogue (S5.6)", fontsize=10.5)

    # 摘要条
    axs = fig.add_subplot(gs[1, 2]); axs.axis("off")
    summary = (
        f"Code 81 audit summary (seed {SEED})\n\n"
        f"H(D) = {H_D:.3f} bits (7 doses, uniform)\n"
        f"I(D; fate) = {I_chain[-1]:.3f} bits end-to-end\n"
        f"Losses L0..L4 = {[round(float(l),2) for l in losses]}\n"
        f"Weakest layer: {weakest_name} ({losses[weakest]:.2f} bits)\n\n"
        f"Gain sweep: I(count) {gain_sweep[0]['I_count']:.2f} -> "
        f"{gain_sweep[-1]['I_count']:.2f} bits (flat)\n"
        f"CV(A) {gain_sweep[0]['cv_A']:.3f} -> {gain_sweep[-1]['cv_A']:.3f} "
        f"(archived 0.002 -> 0.639)\n"
        f"CV(T) pinned <= {max(g['cv_T'] for g in gain_sweep):.3f} (archived 0.016)\n\n"
        f"Amplitude share @ $\\sigma_k$=0.6: {amp_share['sigma_k=0.6']['share']*100:.1f}%\n"
        f"Energy: L*=3.5, p53 L=4 -> counting wins\n\n"
        f"Registered gap: L1 (ATM sensing, no data)"
    )
    axs.text(0.02, 0.98, summary, va="top", fontsize=8.6, family="monospace",
             transform=axs.transAxes)
    fig.suptitle("Code 81: full-pathway p53 cascade audit (L0 dose -> L4 fate)", fontsize=13)
    fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
    log("图已写出:", OUT_PNG)
    log("[代码81] 完成。")


if __name__ == "__main__":
    main()
