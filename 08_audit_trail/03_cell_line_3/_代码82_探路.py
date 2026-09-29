# -*- coding: utf-8 -*-
"""代码82 探路：确认 Mönke 模型振荡区间、轨迹尺度、脉冲检测阈值。用完即弃。"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks, peak_widths

P2025 = dict(A=30.5, P=22, C=1.4, g=2.5, dAM=20, Tm=1.2, TM=4, Tw=1.2, TW=1,
             dA=0.16, dP=0.1, dm=1, dM=2, dw=1.3, dW=2.3,
             kA=0.5, kWA=0.14, kMP=0.15, kPm=1, kPw=1, R=2, Smax=0.2, gamma=9)
P2017 = dict(P2025, Tm=1.0, Tw=1.0)

def rhs(t, y, par, dsb_of_t, s_type):
    ATM, P53, mdm2, Mdm2, wip1, Wip1 = np.maximum(y, 0.0)
    DSB = dsb_of_t(t)
    if s_type == "sat":
        S = par["Smax"] * DSB / (par["gamma"] + DSB)
    else:
        S = par["Smax"] * np.log(DSB / par["gamma"] + 1.0)
    dATM = par["A"] * (ATM**2 / (1 + ATM**2 / par["kA"]**2)) * (1 / (1 + Wip1 / par["kWA"])) \
        - par["dA"] * ATM - par["P"] * ATM * Wip1 + S
    dP53 = par["C"] - par["dP"] * P53 - par["g"] * Mdm2 * P53 / (par["kMP"] + P53) / (1 + par["R"] * ATM)
    dmdm2 = par["Tm"] * P53 / (par["kPm"] + P53) - par["dm"] * mdm2
    dMdm2 = par["TM"] * mdm2 - par["dM"] * Mdm2 - par["dAM"] * ATM * Mdm2
    dwip1 = par["Tw"] * P53 / (par["kPw"] + P53) - par["dw"] * wip1
    dWip1 = par["TW"] * wip1 - par["dW"] * Wip1
    return [dATM, dP53, dmdm2, dMdm2, dwip1, dWip1]

def sim(par, dsb_of_t, s_type, T, y0, dt=0.02):
    t_eval = np.arange(0, T + dt, dt)
    sol = solve_ivp(rhs, (0, T), y0, args=(par, dsb_of_t, s_type),
                    method="LSODA", t_eval=t_eval, rtol=1e-7, atol=1e-10, max_step=0.1)
    return sol.t, sol.y

# 基态（DSB=2 恒定，200h 预平衡）
par = P2025
t, y = sim(par, lambda t: 2.0, "sat", 200.0, [0, 0, 0, 0, 0, 0], dt=0.5)
y0 = y[:, -1]
print("ground state @DSB=2:", np.round(y0, 4))

# R1: 恒定高 DSB 扫描找振荡
for par_name, pp in [("2025", P2025), ("2017", P2017)]:
    for D in [5, 10, 20, 40, 80]:
        t, y = sim(pp, lambda t, D=D: float(D), "sat", 96.0, y0)
        p53 = y[1]
        m = t > 36
        tt, pp53 = t[m], p53[m]
        p10, p99 = np.percentile(pp53, [10, 99])
        pk, props = find_peaks(pp53, prominence=0.2 * (p99 - p10 + 1e-9),
                               height=p10 + 0.25 * (p99 - p10), distance=int(2 / 0.02))
        ipis = np.diff(tt[pk])
        amp = pp53[pk] if len(pk) else []
        print(f"{par_name} DSB={D:>3}: n_pulse={len(pk)} "
              f"T={np.mean(ipis):.2f}+-{np.std(ipis):.2f} " if len(pk) >= 2 else f"{par_name} DSB={D:>3}: n_pulse={len(pk)} ",
              f"p53 range [{pp53.min():.3f},{pp53.max():.3f}]",
              f"amp={np.mean(amp):.3f}" if len(pk) else "")

# R4 探路：log vs sat 高剂量行为
print("--- R4 probe ---")
for s_type in ["sat", "log"]:
    for D in [1, 10, 50, 200, 500]:
        t, y = sim(P2025, lambda t, D=D: float(D), s_type, 96.0, y0)
        p53 = y[1][t > 48]
        rel = (np.percentile(p53, 95) - np.percentile(p53, 5)) / np.median(p53)
        print(f"{s_type} DSB={D:>4}: mean={p53.mean():.3f} rel_osc={rel:.3f}")
