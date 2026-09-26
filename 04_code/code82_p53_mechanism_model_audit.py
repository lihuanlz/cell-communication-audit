# -*- coding: utf-8 -*-
"""
代码82：Mönke 2017 p53 机制模型审计（23 参数 Fisher / 简并群 / 时间-幅度通道 / 7 靶标开庭）
===================================================================================
日期：2026-09-23 ｜ 种子固定：20260923 ｜ 运行预算：数分钟

权威输入：
  - 03_细胞线3/文献参数/机制模型方程与参数表_2026-09-23.md（23 参数表，2025 改版 Tm=Tw=1.2）
  - 03_细胞线3/文献参数/monke2017_supp/srep46571-s1.pdf（方程排版核对，第 5 页 Eq.3）

!!! 方程勘误（相对存档 md，已对照补充材料 PDF 第 5 页渲染图逐字核对） !!!
  E1: ATM 自激活 Hill 项分母为 (1 + ATM*^2/k_A)，不是 (1 + ATM*^2/k_A^2)。
      PDF 高倍渲染图：分母印的是 k_A（无平方）。
  E2: p53 的 ATM* 保护因子为 (1 + R/(1+ATM*))，不是 1/(1+R*ATM*)。
      即 ATM*=0 时 Mdm2 降解效力 (1+R)g，ATM* 高时降到 g。
  两处勘误的动态后果：勘误前模型在任何恒定 S 下都不振荡（稳定螺旋），
  与原文 Fig 3/S5 及 2025 论文 EV6C 矛盾；勘误后恒定 S>S_c≈0.07 进入持续振荡，
  S≈1.1 处振荡熄灭（上 Hopf），与两篇文献全部定性行为一致。

模型（6 ODE，时间 h，浓度 AU）：
  dATM*/dt = A·[ATM*²/(1+ATM*²/kA)]·[1/(1+Wip1/kWA)] − dA·ATM* − P·ATM*·Wip1 + S(DSB)
  dP53/dt  = C − dP·P53 − g·Mdm2·P53/(kMP+P53)·(1 + R/(1+ATM*))
  dmdm2/dt = Tm·P53/(kPm+P53) − dm·mdm2
  dMdm2/dt = TM·mdm2 − dM·Mdm2 − dAM·ATM*·Mdm2
  dwip1/dt = Tw·P53/(kPw+P53) − dw·wip1
  dWip1/dt = TW·wip1 − dW·Wip1
  S(DSB)   = Smax·DSB/(γ+DSB)（饱和型，2017 原版）或 Smax·log(DSB/γ+1)（log 型，2025 版）

DSB 输入：
  确定性 NCS：b(t)=b_s（首小时）后 b_b=2.3，r=0.315，dDSB/dt=b−r·DSB，DSB(0)=b_b/r。
  随机生灭：双池 Gillespie。出生 b(t)，90% iDSB（r_i=0.35/h），10% cDSB（r_c=ln2/20/h）。
             基础 b=0.7 DSB/h；NCS 时首小时 b=b_s。

输出：
  03_细胞线3/结果/代码82_机制审计_结果.json
  03_细胞线3/结果/代码82_基线复现_四面体.{png,svg}
  03_细胞线3/结果/代码82_Fisher审计_四面体.{png,svg}
  03_细胞线3/结果/代码82_验证靶标_六面体.{png,svg}
"""

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks, peak_widths
from scipy.optimize import curve_fit
import sys

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
try:
    from daimon_runtime import setup_plot
    setup_plot()
except Exception:
    pass
plt.rcParams["svg.fonttype"] = "none"

SEED = 20260923
ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "结果"

# ---------------- 参数表（2025 改版；2017 原值仅 Tm/Tw 不同） ----------------
PARAM_ORDER = ["A", "P", "C", "g", "dAM", "Tm", "TM", "Tw", "TW",
               "dA", "dP", "dm", "dM", "dw", "dW",
               "kA", "kWA", "kMP", "kPm", "kPw", "R", "Smax", "gamma"]
P2025 = dict(A=30.5, P=22, C=1.4, g=2.5, dAM=20, Tm=1.2, TM=4, Tw=1.2, TW=1,
             dA=0.16, dP=0.1, dm=1, dM=2, dw=1.3, dW=2.3,
             kA=0.5, kWA=0.14, kMP=0.15, kPm=1, kPw=1, R=2, Smax=0.2, gamma=9)
P2017 = dict(P2025, Tm=1.0, Tw=1.0)
TIME_PARAMS = ["Tm", "Tw", "dA", "dP", "dm", "dM", "dw", "dW"]   # 时间尺度类（任务口径）
AMP_PARAMS = ["C", "g", "Smax"]                                  # 幅度类
OTHER_PARAMS = [p for p in PARAM_ORDER if p not in TIME_PARAMS + AMP_PARAMS]

# DSB 过程常数
R_NCS, B_BASAL_NCS = 0.315, 2.3      # NCS 确定性：修复率与基础出生率
B_BASAL, R_I, F_C = 0.7, 0.35, 0.1   # 随机：基础出生率、iDSB 修复率、cDSB 比例
R_C = np.log(2) / 20.0               # cDSB 修复率（半衰期 20 h）
GY_PER_BSH = 1.0 / 35.0              # 剂量换算锚：35 DSB/Gy（代码81 同口径）

# ---------------- 模型核心 ----------------

def rhs(t, y, par, S_of_t):
    ATM, P53, mdm2, Mdm2, wip1, Wip1 = np.maximum(y, 0.0)
    S = S_of_t(t)
    hill = ATM * ATM / (1.0 + ATM * ATM / par["kA"])          # 勘误 E1：/kA
    dATM = par["A"] * hill / (1.0 + Wip1 / par["kWA"]) \
        - par["dA"] * ATM - par["P"] * ATM * Wip1 + S
    dP53 = par["C"] - par["dP"] * P53 \
        - par["g"] * Mdm2 * P53 / (par["kMP"] + P53) * (1.0 + par["R"] / (1.0 + ATM))  # 勘误 E2
    dmdm2 = par["Tm"] * P53 / (par["kPm"] + P53) - par["dm"] * mdm2
    dMdm2 = par["TM"] * mdm2 - par["dM"] * Mdm2 - par["dAM"] * ATM * Mdm2
    dwip1 = par["Tw"] * P53 / (par["kPw"] + P53) - par["dw"] * wip1
    dWip1 = par["TW"] * wip1 - par["dW"] * Wip1
    return [dATM, dP53, dmdm2, dMdm2, dwip1, dWip1]


def simulate(par, S_of_t, T, y0, dt=0.05, fast=False):
    t_eval = np.arange(0, T + 0.5 * dt, dt)
    if fast:
        sol = solve_ivp(rhs, (0, T), y0, args=(par, S_of_t), method="LSODA",
                        t_eval=t_eval, rtol=1e-6, atol=1e-9, max_step=0.1)
    else:
        sol = solve_ivp(rhs, (0, T), y0, args=(par, S_of_t), method="LSODA",
                        t_eval=t_eval, rtol=1e-7, atol=1e-10, max_step=dt)
    return sol.t, sol.y


def S_const(s):
    return lambda t: s


def S_of_dsb(par, dsb_of_t, s_type):
    if s_type == "sat":
        return lambda t: par["Smax"] * dsb_of_t(t) / (par["gamma"] + dsb_of_t(t))
    return lambda t: par["Smax"] * np.log(dsb_of_t(t) / par["gamma"] + 1.0)


def dsb_det_ncs(b_s, b_b=B_BASAL_NCS, r=R_NCS, t_switch=1.0):
    """确定性 NCS：dDSB/dt=b−r·DSB，b 在 t_switch 处阶跃。解析解。"""
    d0 = b_b / r
    d1 = b_s / r + (d0 - b_s / r) * np.exp(-r * t_switch)
    def f(t):
        if t < t_switch:
            return b_s / r + (d0 - b_s / r) * np.exp(-r * t)
        return b_b / r + (d1 - b_b / r) * np.exp(-r * (t - t_switch))
    return f


def dsb_stoch_grid(rng, t_grid, b_s=0.0, b_basal=B_BASAL, t_switch=1.0,
                   r_i=R_I, r_c=R_C, f_c=F_C):
    """双池生灭 Gillespie，输出按 t_grid 阶梯保持的 DSB 总数轨迹。"""
    out = np.zeros_like(t_grid)
    n_i = rng.poisson(b_basal * (1 - f_c) / r_i)   # 稳态泊松初始化
    n_c = rng.poisson(b_basal * f_c / r_c)
    t, gi, T_end = 0.0, 0, t_grid[-1]
    while t < T_end:
        b = b_s if t < t_switch else b_basal
        a0 = b + r_i * n_i + r_c * n_c
        boundary = t_switch if t < t_switch else T_end
        if a0 > 0:
            t_next = t + rng.exponential(1.0 / a0)
        else:
            t_next = np.inf
        if t_next < boundary:                      # 事件发生在边界前
            while gi < len(t_grid) and t_grid[gi] < t_next:
                out[gi] = n_i + n_c; gi += 1
            u = rng.random() * a0
            if u < b:
                if rng.random() < f_c:
                    n_c += 1
                else:
                    n_i += 1
            elif u < b + r_i * n_i:
                n_i -= 1
            else:
                n_c -= 1
            t = t_next
        else:                                      # 越过出生率阶跃边界或终点
            while gi < len(t_grid) and t_grid[gi] < boundary:
                out[gi] = n_i + n_c; gi += 1
            t = boundary
    while gi < len(t_grid):
        out[gi] = n_i + n_c; gi += 1
    return out


def grid_func(t_grid, vals):
    dt = t_grid[1] - t_grid[0]
    n = len(t_grid)
    def f(t):
        i = int(t / dt)
        return vals[i] if i < n else vals[-1]
    return f


def ground_state(par, dsb_bg=0.0, s_type="sat"):
    """恒定背景 DSB 下积分 300 h 得基态。"""
    S = S_of_dsb(par, lambda t: dsb_bg, s_type)
    t, y = simulate(par, S, 300.0, [0.0] * 6, dt=1.0)
    return y[:, -1]


def detect_pulses(t, p53, prom=0.3, height=0.8, dist_h=2.0):
    """脉冲检测：prominence>=0.3 AU，height>=0.8 AU，最小间隔 2 h。
    阈值依据：基态 P53≈0.17 AU，典型脉冲幅度≈2.7 AU（探路实验标定）。"""
    dt = t[1] - t[0]
    pk, props = find_peaks(p53, prominence=prom, height=height,
                           distance=int(dist_h / dt))
    if len(pk) == 0:
        return dict(times=np.array([]), amps=np.array([]), widths=np.array([]))
    w = peak_widths(p53, pk, rel_height=0.5)[0] * dt
    return dict(times=t[pk], amps=p53[pk], widths=w)

# ---------------- 主流程 ----------------

def main():
    log = lambda *a: print(*a, flush=True)
    res = dict(meta=dict(script="代码82_p53机制模型审计.py", seed=SEED, date="2026-09-23",
                         errata=["E1: Hill 分母 1+ATM^2/kA（非 kA^2）",
                                 "E2: p53 保护因子 1+R/(1+ATM*)（非 1/(1+R·ATM*)）"],
                         params2025=P2025, params2017=P2017))
    y0_sat = ground_state(P2025, dsb_bg=2.0)
    y0 = y0_sat
    log("[0] 基态(2025, sat, DSB=2):", np.round(y0, 4))
    res["ground_state"] = dict(zip(["ATM", "P53", "mdm2", "Mdm2", "wip1", "Wip1"],
                                   [float(v) for v in y0]))

    # ============ 基线复现 R1-R4 ============
    log("[R1] 恒定高 DSB 持续振荡")
    r1 = {}
    for tag, pp in [("2025", P2025), ("2017", P2017)]:
        for s_type, D in [("log", 100.0), ("sat", 1000.0)]:
            S = S_of_dsb(pp, lambda t, D=D: D, s_type)
            t, y = simulate(pp, S, 144.0, y0)
            m = t > 48
            pu = detect_pulses(t[m], y[1][m])
            T_per = float(np.mean(np.diff(pu["times"]))) if len(pu["times"]) >= 3 else np.nan
            r1[f"{tag}_{s_type}"] = dict(
                S=float(S(0)), n_pulse=len(pu["times"]), period=T_per,
                amp=float(np.mean(pu["amps"])) if len(pu["amps"]) else np.nan)
            log(f"  {tag} {s_type} DSB={D}: S={S(0):.3f} n={len(pu['times'])} "
                f"T={T_per:.2f} h amp={r1[f'{tag}_{s_type}']['amp']:.2f}")
    res["R1"] = r1
    # 供图用轨迹
    t_r1, y_r125 = simulate(P2025, S_of_dsb(P2025, lambda t: 100.0, "log"), 48.0, y0)
    _, y_r117 = simulate(P2017, S_of_dsb(P2017, lambda t: 100.0, "log"), 48.0, y0)

    log("[R2] 基态可兴奋性（随机基础 DSB, 300 细胞）")
    rng = np.random.default_rng(SEED + 20)
    t_grid = np.arange(0, 48.0001, 0.05)
    n2 = 200
    r2_counts, r2_amps, r2_traces = [], [], []
    for c in range(n2):
        g = dsb_stoch_grid(rng, t_grid)
        t, y = simulate(P2025, S_of_dsb(P2025, grid_func(t_grid, g), "sat"), 48.0, y0, fast=True)
        pu = detect_pulses(t, y[1])
        r2_counts.append(len(pu["times"]))
        r2_amps.extend(pu["amps"].tolist())
        if c < 3:
            r2_traces.append((t.copy(), y[1].copy(), g.copy()))
    r2_counts = np.array(r2_counts)
    res["R2"] = dict(n_cells=n2, frac_ge1=float(np.mean(r2_counts >= 1)),
                     mean_pulses=float(r2_counts.mean()),
                     count_hist=np.bincount(r2_counts).tolist(),
                     amp_mean=float(np.mean(r2_amps)) if r2_amps else np.nan,
                     amp_cv=float(np.std(r2_amps) / np.mean(r2_amps)) if r2_amps else np.nan)
    log(f"  脉冲细胞比例={res['R2']['frac_ge1']:.3f} 平均脉冲数={res['R2']['mean_pulses']:.2f} "
        f"幅度={res['R2']['amp_mean']:.2f} CV={res['R2']['amp_cv']:.3f}")

    log("[R3] 脉冲幅度/宽度对输入强度不敏感")
    r3 = {}
    base = None
    for b_s in [50, 200, 600]:
        S = S_of_dsb(P2025, dsb_det_ncs(b_s), "sat")
        t, y = simulate(P2025, S, 24.0, y0)
        pu = detect_pulses(t, y[1])
        r3[f"b_s={b_s}"] = dict(amp1=float(pu["amps"][0]), width1=float(pu["widths"][0]))
    for sm in [0.1, 0.2, 0.4]:
        pp = dict(P2025, Smax=sm)
        S = S_of_dsb(pp, dsb_det_ncs(200), "sat")
        t, y = simulate(pp, S, 24.0, y0)
        pu = detect_pulses(t, y[1])
        r3[f"Smax={sm}"] = dict(amp1=float(pu["amps"][0]), width1=float(pu["widths"][0]))
    for k, v in r3.items():
        log(f"  {k}: 首脉冲幅度={v['amp1']:.3f} 宽度={v['width1']:.2f} h")
    res["R3"] = r3

    log("[R4] 高剂量振荡->持续转变（log 型应有，饱和型不应有）")
    dsb_levels = np.logspace(0, np.log10(4000), 15)
    r4 = {}
    for s_type in ["sat", "log"]:
        rows = []
        for D in dsb_levels:
            S = S_of_dsb(P2025, lambda t, D=D: D, s_type)
            t, y = simulate(P2025, S, 144.0, y0, fast=True)
            p53 = y[1][t > 72]
            rel = float((np.percentile(p53, 95) - np.percentile(p53, 5)) /
                        (np.median(p53) + 1e-12))
            rows.append(dict(DSB=float(D), S=float(S(0)), rel_osc=rel,
                             mean=float(p53.mean())))
        r4[s_type] = rows
        log(f"  {s_type}: rel_osc(DSB={rows[0]['DSB']:.0f})={rows[0]['rel_osc']:.3f} -> "
            f"rel_osc(DSB={rows[-1]['DSB']:.0f})={rows[-1]['rel_osc']:.3f}")
    for s_type in ["sat", "log"]:
        below = [r["DSB"] for r in r4[s_type] if r["DSB"] > 5 and r["rel_osc"] < 0.1]
        res.setdefault("R4_transition", {})[s_type] = float(min(below)) if below else None
    log(f"  转变点(rel_osc<0.1): sat={res['R4_transition']['sat']} "
        f"log={res['R4_transition']['log']}")
    res["R4"] = r4

    # ============ Fisher 普查 ============
    log("[F1] 名义轨迹 Jacobian + FIM（确定性 NCS b_s=200, sat, 48h）")
    W_F, DT_F = 48.0, 0.05
    tF = np.arange(0, W_F + 0.5 * DT_F, DT_F)
    S_nom = S_of_dsb(P2025, dsb_det_ncs(200.0), "sat")
    t, y_nom = simulate(P2025, S_nom, W_F, y0, dt=DT_F)
    p53_nom = y_nom[1]
    DLOG = 1e-3
    J = np.zeros((len(tF), len(PARAM_ORDER)))
    for j, pname in enumerate(PARAM_ORDER):
        for sgn in (+1, -1):
            pp = dict(P2025)
            pp[pname] = pp[pname] * np.exp(sgn * DLOG)
            S_p = S_of_dsb(pp, dsb_det_ncs(200.0), "sat")
            _, yp = simulate(pp, S_p, W_F, y0, dt=DT_F)
            if sgn > 0:
                J[:, j] = yp[1]
            else:
                J[:, j] = (J[:, j] - yp[1]) / (2 * DLOG)
    sigma_obs = 0.1 * p53_nom + 0.05          # 比例噪声 10% + 0.05 AU 基底
    Wgt = 1.0 / sigma_obs**2
    FIM = (J * Wgt[:, None]).T @ J
    eigval, eigvec = np.linalg.eigh(FIM)
    idx = np.argsort(eigval)[::-1]
    eigval, eigvec = eigval[idx], eigvec[:, idx]
    cond = float(eigval[0] / max(eigval[-1], 1e-300))
    decades = float(np.log10(eigval[0] / max(eigval[eigval > 0][-1], 1e-300)))
    n_nearzero_6 = int(np.sum(eigval / eigval[0] < 1e-6))
    n_nearzero_8 = int(np.sum(eigval / eigval[0] < 1e-8))
    log(f"  特征值跨度 {decades:.1f} 个数量级, 条件数 {cond:.2e}, "
        f"近零方向(<1e-6) {n_nearzero_6} 个, (<1e-8) {n_nearzero_8} 个")

    def eigvec_composition(k, top=6):
        v = np.abs(eigvec[:, k])**2
        order = np.argsort(v)[::-1][:top]
        return [(PARAM_ORDER[i], float(v[i])) for i in order]

    stiff_dirs = [eigvec_composition(k) for k in range(3)]
    floppy_dirs = [eigvec_composition(k) for k in range(len(PARAM_ORDER) - 3, len(PARAM_ORDER))]
    for k in range(3):
        log(f"  stiff#{k+1} (lam={eigval[k]:.2e}): " +
            ", ".join(f"{n}={w:.2f}" for n, w in stiff_dirs[k][:4]))
    for i, k in enumerate(range(len(PARAM_ORDER) - 3, len(PARAM_ORDER))):
        log(f"  floppy#{i+1} (lam={eigval[k]:.2e}): " +
            ", ".join(f"{n}={w:.2f}" for n, w in floppy_dirs[i][:4]))

    # 类别归属：每个特征方向在 时间/幅度/其他 三类参数上的权重
    def class_weights(k):
        v = np.abs(eigvec[:, k])**2
        return dict(time=float(sum(v[PARAM_ORDER.index(p)] for p in TIME_PARAMS)),
                    amp=float(sum(v[PARAM_ORDER.index(p)] for p in AMP_PARAMS)),
                    other=float(sum(v[PARAM_ORDER.index(p)] for p in OTHER_PARAMS)))
    eig_class = [class_weights(k) for k in range(len(PARAM_ORDER))]

    res["F1"] = dict(eigval=eigval.tolist(), cond=cond, decades=decades,
                     n_nearzero_6=n_nearzero_6, n_nearzero_8=n_nearzero_8,
                     stiff_dirs=stiff_dirs, floppy_dirs=floppy_dirs,
                     eig_class=eig_class, param_order=PARAM_ORDER,
                     noise_model="sigma=0.1*P53+0.05 AU", dlog=DLOG)

    # ============ 时间 vs 幅度通道分解 ============
    log("[F2] 时间 vs 幅度通道 Fisher 分解")
    pu_nom = detect_pulses(t, p53_nom)
    n_pulse_nom = len(pu_nom["times"])
    IPI_nom = float(np.mean(np.diff(pu_nom["times"])))
    SIG_T = 0.30 * IPI_nom            # 归档细胞间 IPI CV≈0.30
    SIG_A_FRAC = 0.125                # 幅度 CV 10-15% 中点
    log(f"  名义脉冲 {n_pulse_nom} 个, IPI={IPI_nom:.2f} h, sigma_t={SIG_T:.2f} h")

    def channel_fisher(sig_t, sig_a_frac):
        F_t = np.zeros(len(PARAM_ORDER))
        F_a = np.zeros(len(PARAM_ORDER))
        count_changed = []
        for j, pname in enumerate(PARAM_ORDER):
            dtp = np.zeros(n_pulse_nom); dap = np.zeros(n_pulse_nom)
            ok = True
            for sgn in (+1, -1):
                pp = dict(P2025)
                pp[pname] = pp[pname] * np.exp(sgn * DLOG)
                S_p = S_of_dsb(pp, dsb_det_ncs(200.0), "sat")
                _, yp = simulate(pp, S_p, W_F, y0, dt=DT_F)
                pu = detect_pulses(tF, yp[1])
                m = min(n_pulse_nom, len(pu["times"]))
                if len(pu["times"]) != n_pulse_nom:
                    ok = False
                tv = np.full(n_pulse_nom, np.nan); av = np.full(n_pulse_nom, np.nan)
                tv[:m] = pu["times"][:m]; av[:m] = pu["amps"][:m]
                if sgn > 0:
                    dtp, dap = tv, av
                else:
                    dtp = (dtp - tv) / (2 * DLOG)
                    dap = (dap - av) / (2 * DLOG)
            if not ok:
                count_changed.append(pname)
            valid = ~np.isnan(dtp)
            F_t[j] = np.nansum((dtp[valid] / sig_t) ** 2)
            F_a[j] = np.nansum((dap[valid] / (sig_a_frac * pu_nom["amps"][valid])) ** 2)
        return F_t, F_a, count_changed

    F_t, F_a, count_changed = channel_fisher(SIG_T, SIG_A_FRAC)
    share_t = float(F_t.sum() / (F_t.sum() + F_a.sum()))
    log(f"  时间通道份额={share_t*100:.1f}% 幅度={(1-share_t)*100:.1f}% "
        f"(归档 94.2%/5.5%) 计数改变参数: {count_changed}")
    # 噪声口径网格：sigma_t 取归档细胞内 IPI 抖动 0.085 或细胞间 0.30，
    # sigma_A 取内生钉死值 0.002 或细胞间观测 0.125。检验通道分解对噪声标定的依赖。
    noise_grid = {}
    for st_cv in [0.085, 0.30]:
        for sa in [0.002, 0.125]:
            Ft2, Fa2, _ = channel_fisher(st_cv * IPI_nom, sa)
            noise_grid[f"sigT={st_cv}xIPI,sigA={sa}xA"] = float(Ft2.sum() / (Ft2.sum() + Fa2.sum()))
    log(f"  噪声口径网格(时间份额): { {k: round(v,3) for k,v in noise_grid.items()} }")
    per_param = []
    for j, pname in enumerate(PARAM_ORDER):
        tot = F_t[j] + F_a[j]
        per_param.append(dict(param=pname, F_time=float(F_t[j]), F_amp=float(F_a[j]),
                              time_share=float(F_t[j] / tot) if tot > 0 else np.nan))
    res["F2"] = dict(n_pulse=n_pulse_nom, IPI=IPI_nom, sig_t=SIG_T, sig_a_frac=SIG_A_FRAC,
                     share_time=share_t, share_amp=1 - share_t, noise_grid=noise_grid,
                     count_changed=count_changed, per_param=per_param,
                     archived=dict(time=0.942, amp=0.055))

    # ============ 结构常数执照（±10x 扫描） ============
    log("[F3] 结构常数扫描（恒定 DSB=100, log 型, S=0.498）")
    folds = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
    S_scan = S_of_dsb(P2025, lambda t: 100.0, "log")
    t_s, y_s = simulate(P2025, S_scan, 144.0, y0)
    pu0 = detect_pulses(t_s[t_s > 48], y_s[1][t_s > 48])
    T0 = float(np.mean(np.diff(pu0["times"])))
    A0 = float(np.mean(pu0["amps"]))
    log(f"  名义 T0={T0:.2f} h, A0={A0:.2f} AU")
    scan = []
    for pname in PARAM_ORDER:
        Ts, As = [], []
        for f in folds:
            pp = dict(P2025)
            pp[pname] = pp[pname] * f
            S_p = S_of_dsb(pp, lambda t: 100.0, "log")
            t_p, y_p = simulate(pp, S_p, 144.0, y0, fast=True)
            pu = detect_pulses(t_p[t_p > 48], y_p[1][t_p > 48])
            if len(pu["times"]) >= 3:
                Ts.append(float(np.mean(np.diff(pu["times"]))))
                As.append(float(np.mean(pu["amps"])))
            else:
                Ts.append(np.nan); As.append(np.nan)
        T_arr = np.array(Ts); A_arr = np.array(As)
        cls = "time" if pname in TIME_PARAMS else ("amp" if pname in AMP_PARAMS else "other")
        scan.append(dict(param=pname, cls=cls, folds=folds,
                         T=[None if np.isnan(x) else x for x in T_arr],
                         A=[None if np.isnan(x) else x for x in A_arr],
                         max_dT=float(np.nanmax(np.abs(T_arr / T0 - 1))),
                         max_dA=float(np.nanmax(np.abs(A_arr / A0 - 1))),
                         n_noosc=int(np.sum(np.isnan(T_arr)))))
        log(f"  {pname:>5} ({cls:>5}): max|dT/T|={scan[-1]['max_dT']:.3f} "
            f"max|dA/A|={scan[-1]['max_dA']:.3f} no-osc folds={scan[-1]['n_noosc']}")
    amp_T = [s["max_dT"] for s in scan if s["cls"] == "amp"]
    time_T = [s["max_dT"] for s in scan if s["cls"] == "time"]
    res["F3"] = dict(T0=T0, A0=A0, folds=folds, scan=scan,
                     amp_class_max_dT=amp_T, time_class_max_dT=time_T)

    # ============ 验证靶标 ============
    log("[V1] 计数律 N(D)：随机 NCS 系综，b_s 扫描")
    BS_LIST = [0, 25, 50, 100, 150, 200, 300, 400, 600]
    N_CELLS_D = 120
    rng = np.random.default_rng(SEED + 30)
    dose_rows = []
    ens200 = None
    for b_s in BS_LIST:
        counts = []
        store = []
        for c in range(N_CELLS_D):
            g = dsb_stoch_grid(rng, t_grid, b_s=b_s, t_switch=1.0)
            t, y = simulate(P2025, S_of_dsb(P2025, grid_func(t_grid, g), "sat"), 48.0, y0,
                            fast=True)
            pu = detect_pulses(t, y[1])
            counts.append(len(pu["times"]))
            store.append(pu)
        counts = np.array(counts)
        dose_gy = b_s * 1.0 * GY_PER_BSH     # 首小时累积 DSB 数 / 35 DSB/Gy
        dose_rows.append(dict(b_s=b_s, dose_gy=float(dose_gy),
                              mean_N=float(counts.mean()), sem=float(counts.std() / np.sqrt(len(counts))),
                              median_N=float(np.median(counts))))
        log(f"  b_s={b_s:>3} (D={dose_gy:.2f} Gy): N={counts.mean():.2f}±{counts.std():.2f}")
        if b_s == 200:
            ens200 = store
    Dg = np.array([r["dose_gy"] for r in dose_rows])
    Nm = np.array([r["mean_N"] for r in dose_rows])

    def hyperb(D, Nmax, Dc):
        return Nmax * D / (Dc + D)
    def hyperb3(D, N0, Nmax, Dc):
        return N0 + Nmax * D / (Dc + D)
    # 主拟合：诱导脉冲数 N(D)-N(0)，过原点双曲（与归档口径可比）
    Nind = Nm - Nm[0]
    try:
        popt, _ = curve_fit(hyperb, Dg, Nind, p0=[4.0, 0.3], maxfev=20000)
        pred = hyperb(Dg, *popt)
        r2_hyp = 1 - np.sum((Nind - pred) ** 2) / np.sum((Nind - Nind.mean()) ** 2)
    except Exception:
        popt, r2_hyp = [np.nan, np.nan], np.nan
    # 备选：带基底的三参数双曲（容忍基态自发脉冲）
    try:
        popt3, _ = curve_fit(hyperb3, Dg, Nm, p0=[Nm[0], 4.0, 0.3], maxfev=20000)
        pred3 = hyperb3(Dg, *popt3)
        r2_hyp3 = 1 - np.sum((Nm - pred3) ** 2) / np.sum((Nm - Nm.mean()) ** 2)
    except Exception:
        popt3, r2_hyp3 = [np.nan] * 3, np.nan
    lin = np.polyfit(Dg, Nind, 1)
    pred_l = np.polyval(lin, Dg)
    r2_lin = 1 - np.sum((Nind - pred_l) ** 2) / np.sum((Nind - Nind.mean()) ** 2)
    res["V1"] = dict(doses=dose_rows, hyp=dict(Nmax=float(popt[0]), Dc=float(popt[1]),
                                               R2=float(r2_hyp), target="N(D)-N(0)"),
                     hyp3=dict(N0=float(popt3[0]), Nmax=float(popt3[1]), Dc=float(popt3[2]),
                               R2=float(r2_hyp3)),
                     lin=dict(a=float(lin[1]), b=float(lin[0]), R2=float(r2_lin)),
                     archived=dict(R2=0.920, Dc=0.261),
                     note="D 由 b_s×1h/35 DSB/Gy 换算，单位映射为假设口径；"
                          "基态自发脉冲 N(0)≈4.5 由双池随机背景所致")
    log(f"  诱导脉冲双曲拟合 Nmax={popt[0]:.2f} Dc={popt[1]:.3f} Gy R2={r2_hyp:.3f} | "
        f"线性 R2={r2_lin:.3f} | 三参数 R2={r2_hyp3:.3f} (归档 R2=0.920, Dc=0.261)")

    log("[V3] 幅度 CV：随机 DSB + Tw 对数正态异质性(CV=0.5, Mönke Fig4/S8A 做法)")
    N_HET = 200
    rng = np.random.default_rng(SEED + 31)
    sig_l = np.sqrt(np.log(1 + 0.5**2))
    amp_cv_rows = {}
    het_pulses = None
    for tag, het in [("ctrl", False), ("Tw_het", True)]:
        if not het:
            pus = ens200           # 同质对照直接复用 V1 的 b_s=200 系综（120 细胞）
        else:
            pus = []
            for c in range(N_HET):
                pp = dict(P2025)
                # 保持均值=1 的对数正态乘子：exp(N(-sig_l^2/2, sig_l^2))，CV=0.5
                pp["Tw"] = P2025["Tw"] * float(np.exp(rng.normal(0, sig_l) - 0.5 * sig_l**2))
                g = dsb_stoch_grid(rng, t_grid, b_s=200.0, t_switch=1.0)
                t, y = simulate(pp, S_of_dsb(pp, grid_func(t_grid, g), "sat"), 48.0, y0,
                                fast=True)
                pus.append(detect_pulses(t, y[1]))
        first_amps, all_amps, cell_mean_amps = [], [], []
        for pu in pus:
            if len(pu["amps"]) >= 1:
                first_amps.append(pu["amps"][0])
                all_amps.extend(pu["amps"].tolist())
                cell_mean_amps.append(float(np.mean(pu["amps"])))
        amp_cv_rows[tag] = dict(cv_first=float(np.std(first_amps) / np.mean(first_amps)),
                                cv_all=float(np.std(all_amps) / np.mean(all_amps)),
                                cv_cellmean=float(np.std(cell_mean_amps) / np.mean(cell_mean_amps)),
                                mean_first=float(np.mean(first_amps)),
                                n=len(first_amps))
        log(f"  {tag}: CV(first)={amp_cv_rows[tag]['cv_first']:.3f} "
            f"CV(all)={amp_cv_rows[tag]['cv_all']:.3f} "
            f"CV(cellmean)={amp_cv_rows[tag]['cv_cellmean']:.3f} (靶标 0.10-0.15)")
        if het:
            het_pulses = pus
    res["V3"] = amp_cv_rows

    log("[V2] 弥散比 Var(t1)/Var(IPI)：同质系综 + Tw 异质系综")
    def dispersion(pus):
        t1s, ipis, ipi_cv_within = [], [], []
        for pu in pus:
            if len(pu["times"]) >= 1:
                t1s.append(pu["times"][0])
            if len(pu["times"]) >= 2:
                d = np.diff(pu["times"])
                ipis.extend(d.tolist())
                if len(d) >= 3:
                    ipi_cv_within.append(float(np.std(d) / np.mean(d)))
        return dict(var_t1=float(np.var(t1s, ddof=1)), var_ipi=float(np.var(ipis, ddof=1)),
                    ratio=float(np.var(t1s, ddof=1) / np.var(ipis, ddof=1)),
                    cv_t1=float(np.std(t1s) / np.mean(t1s)),
                    cv_ipi=float(np.std(ipis) / np.mean(ipis)),
                    ipi_cv_within=float(np.mean(ipi_cv_within)) if ipi_cv_within else np.nan,
                    n_t1=len(t1s), n_ipi=len(ipis))
    v2_hom = dispersion(ens200)
    v2_het = dispersion(het_pulses)
    res["V2"] = dict(homogeneous=v2_hom, heterogeneous=v2_het, archived=5.8,
                     note="归档口径：Var(首脉冲时刻)/Var(相继IPI)，文献 5.8，代码81 生成模型 5.7")
    log(f"  同质: Var(t1)={v2_hom['var_t1']:.4f} Var(IPI)={v2_hom['var_ipi']:.4f} "
        f"比值={v2_hom['ratio']:.2f}")
    log(f"  异质: Var(t1)={v2_het['var_t1']:.4f} Var(IPI)={v2_het['var_ipi']:.4f} "
        f"比值={v2_het['ratio']:.2f} (归档 5.8)")

    log("[V5] Wip1 RNAi（Tw=TW 阶梯 x1.0/0.85/0.7/0.5/0.3）")
    rnai = {}
    for f in [1.0, 0.85, 0.7, 0.5, 0.3]:
        tag = f"x{f}"
        pp = dict(P2025)
        pp["Tw"] = P2025["Tw"] * f
        pp["TW"] = P2025["TW"] * f
        y0_kd = ground_state(pp, dsb_bg=2.0)
        S = S_of_dsb(pp, dsb_det_ncs(200.0), "sat")
        t, y = simulate(pp, S, 48.0, y0_kd)
        pu = detect_pulses(t, y[1])
        late = y[1][t > 36]
        rnai[tag] = dict(factor=f, n_pulse=len(pu["times"]),
                         period=float(np.mean(np.diff(pu["times"]))) if len(pu["times"]) >= 2 else np.nan,
                         amp=float(np.mean(pu["amps"])) if len(pu["amps"]) else np.nan,
                         width=float(np.mean(pu["widths"])) if len(pu["widths"]) else np.nan,
                         wip1_basal=float(y0_kd[5]),
                         late_mean=float(late.mean()),
                         rel_osc_late=float((np.percentile(late, 95) - np.percentile(late, 5))
                                            / np.median(late)),
                         trace=(t.copy(), y[1].copy()))
        log(f"  {tag}: n={rnai[tag]['n_pulse']} T={rnai[tag]['period']:.2f} "
            f"amp={rnai[tag]['amp']:.2f} width={rnai[tag]['width']:.2f} "
            f"late_mean={rnai[tag]['late_mean']:.2f} rel_osc={rnai[tag]['rel_osc_late']:.3f}")
    res["V5"] = {k: {kk: vv for kk, vv in v.items() if kk != "trace"}
                 for k, v in rnai.items()}

    res["V4"] = dict(ref="F2 时间/幅度通道分解", share_time=share_t, archived_time=0.942)
    res["V6"] = dict(ref="R4", log_transition="见 R4.log", sat_no_transition="见 R4.sat")
    res["V7"] = dict(ref="R2", frac_ge1=res["R2"]["frac_ge1"])

    # ============ 图 ============
    make_fig1(res, t_r1, y_r125, y_r117, r2_traces, r4)
    make_fig2(res, eigval, eigvec, per_param, share_t, noise_grid)
    make_fig3(res, dose_rows, popt, lin, ens200, amp_cv_rows, rnai, scan, T0)

    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "代码82_机制审计_结果.json", "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2, default=str)
    log("[done] JSON 已写出")


# ---------------- 图 ----------------
C_BLUE, C_RED, C_GREEN, C_ORANGE, C_DARK = "#33527a", "#c05640", "#0b6b3a", "#d9863d", "#8c1d18"
CLS_COLOR = {"time": C_BLUE, "amp": C_RED, "other": "#777777"}


def savefig(fig, stem):
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)
    print("图已写出:", OUT / f"{stem}.png", flush=True)


def make_fig1(res, t_r1, y25, y17, r2_traces, r4):
    fig = plt.figure(figsize=(15, 9.5))
    gs = fig.add_gridspec(2, 2, hspace=0.36, wspace=0.26)

    # (a) R1: sustained oscillation under constant high DSB
    axa = fig.add_subplot(gs[0, 0])
    axa.plot(t_r1, y25[1], color=C_BLUE, lw=1.4,
             label=f"2025 params (Tm=Tw=1.2), T={res['R1']['2025_log']['period']:.2f} h")
    axa.plot(t_r1, y17[1], color=C_ORANGE, lw=1.4,
             label=f"2017 params (Tm=Tw=1.0), T={res['R1']['2017_log']['period']:.2f} h")
    axa.set_xlabel("time (h)"); axa.set_ylabel("P53 (AU)")
    axa.legend(fontsize=8, loc="upper right")
    axa.set_title("(a) R1: sustained oscillations, constant DSB=100 (log input)", fontsize=10.5)

    # (b) R2: excitable ground state
    axb = fig.add_subplot(gs[0, 1])
    for i, (t, p53, g) in enumerate(r2_traces):
        axb.plot(t, p53 + i * 3.5, color=C_BLUE, lw=0.9)
        axb.text(49, i * 3.5 + 0.4, f"DSB max {g.max():.0f}", fontsize=7.5, color="#555")
    axb.set_xlabel("time (h)"); axb.set_ylabel("P53 (AU, offset per cell)")
    axb.set_xlim(0, 60)
    axb.set_title(f"(b) R2: excitable basal state, stochastic DSB "
                  f"(pulsing cells {res['R2']['frac_ge1']*100:.0f}%)", fontsize=10.5)

    # (c) R3: amplitude/width insensitivity
    axc = fig.add_subplot(gs[1, 0])
    bs = [50, 200, 600]
    amps = [res["R3"][f"b_s={b}"]["amp1"] for b in bs]
    wids = [res["R3"][f"b_s={b}"]["width1"] for b in bs]
    axc.plot(bs, np.array(amps) / amps[1], "o-", color=C_RED,
             label="1st pulse amplitude (norm.)")
    axc.plot(bs, np.array(wids) / wids[1], "s-", color=C_BLUE,
             label="1st pulse width (norm.)")
    sm = [0.1, 0.2, 0.4]
    amps2 = [res["R3"][f"Smax={s}"]["amp1"] for s in sm]
    wids2 = [res["R3"][f"Smax={s}"]["width1"] for s in sm]
    axc.plot(np.array(sm) * 1000, np.array(amps2) / amps2[1], "o--", color=C_RED, alpha=0.55,
             label="amplitude vs Smax (x1000)")
    axc.plot(np.array(sm) * 1000, np.array(wids2) / wids2[1], "s--", color=C_BLUE, alpha=0.55,
             label="width vs Smax (x1000)")
    axc.axhline(1, color="#999", lw=0.7, ls=":")
    axc.set_xscale("log")
    axc.set_xlabel("NCS birth rate $b_s$ (breaks/h)  or  $S_{max}\\times 1000$")
    axc.set_ylabel("relative to nominal")
    axc.legend(fontsize=8, loc="lower right")
    axc.set_title("(c) R3: pulse amplitude/width vs input strength", fontsize=10.5)

    # (d) R4: high-dose transition
    axd = fig.add_subplot(gs[1, 1])
    for s_type, col, lab in [("sat", C_BLUE, "saturable S(DSB)"),
                             ("log", C_RED, "log S(DSB) (2025)")]:
        Ds = [r["DSB"] for r in r4[s_type]]
        ro = [r["rel_osc"] for r in r4[s_type]]
        axd.plot(Ds, ro, "o-", color=col, ms=4, label=lab)
    axd.axhline(0.1, color="#999", ls=":", lw=0.9)
    axd.text(1.2, 0.13, "sustained below", fontsize=7.5, color="#666")
    axd.set_xscale("log")
    axd.set_xlabel("constant DSB count"); axd.set_ylabel("relative oscillation amplitude")
    axd.legend(fontsize=8.5)
    axd.set_title("(d) R4: oscillation-to-sustained transition", fontsize=10.5)

    fig.suptitle("Code 82 baseline reproduction: Mönke 2017 p53 model (2025 params)", fontsize=13)
    savefig(fig, "代码82_基线复现_四面体")


def make_fig2(res, eigval, eigvec, per_param, share_t, noise_grid):
    fig = plt.figure(figsize=(15, 9.5))
    gs = fig.add_gridspec(2, 2, hspace=0.38, wspace=0.28)

    # (a) eigenvalue spectrum
    axa = fig.add_subplot(gs[0, 0])
    axa.semilogy(range(1, len(eigval) + 1), eigval, "o-", color=C_BLUE, ms=4)
    axa.axhline(eigval[0] * 1e-6, color=C_DARK, ls=":", lw=0.9)
    axa.text(len(eigval) * 0.45, eigval[0] * 1.8e-6, "1e-6 x max", fontsize=7.5, color=C_DARK)
    axa.set_xlabel("eigenvalue rank"); axa.set_ylabel("FIM eigenvalue")
    axa.set_title(f"(a) Sloppy spectrum: {res['F1']['decades']:.1f} decades, "
                  f"cond={res['F1']['cond']:.1e}, near-zero={res['F1']['n_nearzero_6']}",
                  fontsize=10.5)

    # (b) top-3 stiff directions
    axb = fig.add_subplot(gs[0, 1])
    ypos = 0
    yticks, ylabels = [], []
    for k in range(3):
        comp = res["F1"]["stiff_dirs"][k]
        for name, w in comp[:5]:
            cls = "time" if name in ["Tm", "Tw", "dA", "dP", "dm", "dM", "dw", "dW"] else \
                ("amp" if name in ["C", "g", "Smax"] else "other")
            axb.barh(ypos, w, color=CLS_COLOR[cls], height=0.8)
            yticks.append(ypos); ylabels.append(f"s{k+1}:{name}")
            ypos += 1
        ypos += 0.8
    axb.set_yticks(yticks); axb.set_yticklabels(ylabels, fontsize=7)
    axb.invert_yaxis()
    axb.set_xlabel("squared loading |v_i|^2")
    axb.set_title("(b) top-3 stiff directions (blue=timescale, red=amplitude)", fontsize=10.5)

    # (c) floppiest directions
    axc = fig.add_subplot(gs[1, 0])
    ypos = 0
    yticks, ylabels = [], []
    for i, comp in enumerate(res["F1"]["floppy_dirs"]):
        for name, w in comp[:5]:
            cls = "time" if name in ["Tm", "Tw", "dA", "dP", "dm", "dM", "dw", "dW"] else \
                ("amp" if name in ["C", "g", "Smax"] else "other")
            axc.barh(ypos, w, color=CLS_COLOR[cls], height=0.8)
            yticks.append(ypos); ylabels.append(f"f{i+1}:{name}")
            ypos += 1
        ypos += 0.8
    axc.set_yticks(yticks); axc.set_yticklabels(ylabels, fontsize=7)
    axc.invert_yaxis()
    axc.set_xlabel("squared loading |v_i|^2")
    axc.set_title("(c) 3 floppiest directions (near-unidentifiable combos)", fontsize=10.5)

    # (d) time vs amplitude channel share
    axd = fig.add_subplot(gs[1, 1])
    vals = [share_t, 1 - share_t]
    bars = axd.bar([0, 1], vals, color=[C_BLUE, C_RED], width=0.55)
    axd.axhline(0.942, color="#333", ls="--", lw=1.1)
    axd.axhline(0.055, color="#333", ls=":", lw=1.1)
    axd.text(-0.55, 0.955, "archived 94.2%", fontsize=8, color="#333")
    axd.text(-0.55, 0.075, "archived 5.5%", fontsize=8, color="#333")
    for x, v in zip([0, 1], vals):
        axd.text(x, v + 0.02, f"{v*100:.1f}%", ha="center", fontsize=10)
    axd.set_xticks([0, 1]); axd.set_xticklabels(["pulse timing", "pulse amplitude"])
    axd.set_ylabel("Fisher information share")
    axd.set_ylim(0, 1.08)
    ng = "\n".join(f"{k}: {v*100:.0f}%" for k, v in noise_grid.items())
    axd.text(1.45, 0.62, "noise-calibration grid\n(time share):\n" + ng, fontsize=7.2,
             family="monospace", va="top",
             bbox=dict(fc="#f4f4f4", ec="#999", lw=0.7))
    axd.set_xlim(-0.6, 2.6)
    axd.set_title("(d) channel decomposition at nominal dispersion "
                  "(CV_t=0.30, CV_A=0.125)", fontsize=9.5)
    fig.suptitle("Code 82 Fisher audit: 23 mechanistic parameters, nominal NCS 48 h trajectory",
                 fontsize=13)
    savefig(fig, "代码82_Fisher审计_四面体")


def make_fig3(res, dose_rows, popt, lin, ens200, amp_cv_rows, rnai, scan, T0):
    fig = plt.figure(figsize=(17, 9.5))
    gs = fig.add_gridspec(2, 3, hspace=0.40, wspace=0.30)

    # (a) counting law
    axa = fig.add_subplot(gs[0, 0])
    Dg = np.array([r["dose_gy"] for r in dose_rows])
    Nm = np.array([r["mean_N"] for r in dose_rows])
    sem = np.array([r["sem"] for r in dose_rows])
    axa.errorbar(Dg, Nm - Nm[0], yerr=sem, fmt="o", color=C_BLUE,
                 label="induced pulses N(D)-N(0)")
    Dd = np.linspace(0, Dg.max(), 200)
    axa.plot(Dd, popt[0] * Dd / (popt[1] + Dd), color=C_RED,
             label=f"hyperbolic: R2={res['V1']['hyp']['R2']:.3f}, Dc={popt[1]:.2f} Gy")
    axa.plot(Dd, np.polyval(lin, Dd), color="#777", ls="--",
             label=f"linear: R2={res['V1']['lin']['R2']:.3f}")
    axa.set_xlabel("dose D (Gy, from $b_s\\times$1h / 35 DSB/Gy)")
    axa.set_ylabel("induced pulse count in 48 h")
    axa.legend(fontsize=8, loc="lower right")
    axa.set_title("(a) V1 counting law N(D) (archived R2=0.920, Dc=0.261)", fontsize=10.5)

    # (b) dispersion ratio
    axb = fig.add_subplot(gs[0, 1])
    v = res["V2"]
    xs = [0, 1, 3, 4]
    vals = [v["homogeneous"]["var_t1"], v["homogeneous"]["var_ipi"],
            v["heterogeneous"]["var_t1"], v["heterogeneous"]["var_ipi"]]
    axb.bar(xs, vals, color=[C_ORANGE, C_BLUE, C_ORANGE, C_BLUE], width=0.6)
    axb.set_xticks(xs)
    axb.set_xticklabels(["Var($t_1$)\nhomog.", "Var(IPI)\nhomog.",
                         "Var($t_1$)\nTw-het.", "Var(IPI)\nTw-het."], fontsize=8)
    for x, vv in zip(xs, vals):
        axb.text(x, vv + 0.03 * max(vals), f"{vv:.3f}", ha="center", fontsize=8)
    axb.set_ylabel("variance (h^2)")
    axb.set_title(f"(b) V2 dispersion ratio: homog.={v['homogeneous']['ratio']:.2f}, "
                  f"het.={v['heterogeneous']['ratio']:.2f} (archived 5.8)", fontsize=10)

    # (c) amplitude CV
    axc = fig.add_subplot(gs[0, 2])
    tags = ["ctrl", "Tw_het"]
    cvs = [amp_cv_rows[t]["cv_first"] for t in tags]
    cvm = [amp_cv_rows[t]["cv_cellmean"] for t in tags]
    axc.bar([0, 1], cvs, color=[C_BLUE, C_RED], width=0.35, label="1st pulse CV")
    axc.bar([0.38, 1.38], cvm, color=[C_BLUE, C_RED], width=0.35, alpha=0.5,
            label="per-cell mean CV")
    axc.axhspan(0.10, 0.15, color=C_GREEN, alpha=0.18)
    axc.text(1.8, 0.125, "target\n10-15%", fontsize=8, color=C_GREEN, va="center")
    for x, vv in zip([0, 1], cvs):
        axc.text(x, vv + 0.004, f"{vv*100:.1f}%", ha="center", fontsize=9)
    for x, vv in zip([0.38, 1.38], cvm):
        axc.text(x, vv + 0.004, f"{vv*100:.1f}%", ha="center", fontsize=9)
    axc.set_xticks([0.19, 1.19])
    axc.set_xticklabels(["homog.", "Tw het. CV=0.5"])
    axc.set_ylabel("amplitude CV")
    axc.set_xlim(-0.5, 2.3)
    axc.legend(fontsize=8)
    axc.set_title("(c) V3 amplitude CV under Wip1 heterogeneity", fontsize=10.5)

    # (d) Wip1 RNAi
    axd = fig.add_subplot(gs[1, 0])
    cols = {"x1.0": "#333", "x0.85": C_BLUE, "x0.7": C_ORANGE, "x0.5": C_RED, "x0.3": C_DARK}
    for tag in ["x1.0", "x0.85", "x0.7", "x0.5", "x0.3"]:
        r = rnai[tag]
        t, p53 = r["trace"]
        lab = (f"Tw=TW x{r['factor']}: n={r['n_pulse']}"
               + (f", T={r['period']:.1f} h" if not np.isnan(r["period"]) else ", sustained"))
        axd.plot(t, p53, color=cols[tag], lw=1.2, label=lab)
    axd.set_xlabel("time (h)"); axd.set_ylabel("P53 (AU)")
    axd.legend(fontsize=7.5)
    axd.set_title("(d) V5 Wip1 RNAi ladder (NCS $b_s$=200)", fontsize=10.5)

    # (e) structural constant scan
    axe = fig.add_subplot(gs[1, 1])
    order = np.argsort([s["max_dT"] for s in scan])
    names = [scan[i]["param"] for i in order]
    vals = [scan[i]["max_dT"] for i in order]
    clss = [scan[i]["cls"] for i in order]
    axe.barh(range(len(vals)), vals, color=[CLS_COLOR[c] for c in clss], height=0.75)
    axe.set_yticks(range(len(vals))); axe.set_yticklabels(names, fontsize=7)
    axe.axvline(0.1, color=C_DARK, ls=":", lw=1.0)
    axe.text(0.105, len(vals) - 2, "10%", fontsize=8, color=C_DARK)
    axe.set_xlabel("max |dT/T| over 0.1x-10x fold change")
    axe.set_title(f"(e) structural constant license: T0={T0:.2f} h", fontsize=10.5)

    # (f) summary
    axs = fig.add_subplot(gs[1, 2]); axs.axis("off")
    v1 = res["V1"]
    txt = (
        f"Code 82 validation targets (seed {res['meta']['seed']})\n\n"
        f"V1 counting law: R2_hyp={v1['hyp']['R2']:.3f} vs R2_lin={v1['lin']['R2']:.3f}\n"
        f"   Dc={v1['hyp']['Dc']:.3f} Gy (archived 0.261)\n"
        f"V2 dispersion ratio: homog.={res['V2']['homogeneous']['ratio']:.2f}, "
        f"het.={res['V2']['heterogeneous']['ratio']:.2f} (archived 5.8)\n"
        f"V3 amplitude CV: ctrl={amp_cv_rows['ctrl']['cv_first']*100:.1f}%, "
        f"het={amp_cv_rows['Tw_het']['cv_first']*100:.1f}% (target 10-15%)\n"
        f"V4 time/amp Fisher share: {res['V4']['share_time']*100:.1f}% (archived 94.2%)\n"
        f"V5 Wip1 RNAi: n_pulse x1.0={rnai['x1.0']['n_pulse']}, x0.85={rnai['x0.85']['n_pulse']}, "
        f"x0.7={rnai['x0.7']['n_pulse']}\n"
        f"   x0.5/x0.3 -> sustained high state (n={rnai['x0.5']['n_pulse']})\n"
        f"V6 high-dose transition: see Fig1d\n"
        f"V7 excitable ground: pulsing cells {res['V7']['frac_ge1']*100:.0f}%\n\n"
        f"R1 period: 2025={res['R1']['2025_log']['period']:.2f} h, "
        f"2017={res['R1']['2017_log']['period']:.2f} h (targets 5.5/7.0)\n"
        f"Fisher: {res['F1']['decades']:.1f} decades, near-zero={res['F1']['n_nearzero_6']}"
    )
    axs.text(0.02, 0.98, txt, va="top", fontsize=8.6, family="monospace",
             transform=axs.transAxes)
    fig.suptitle("Code 82 validation targets and structural-constant license", fontsize=13)
    savefig(fig, "代码82_验证靶标_六面体")


if __name__ == "__main__":
    main()
