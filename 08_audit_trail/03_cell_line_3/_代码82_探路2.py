# -*- coding: utf-8 -*-
"""探路2：直接扫恒定 S（绕过 S(DSB)），定位 Hopf 振荡窗口；并测试 NCS 衰减输入下的持续振荡。"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

P2025 = dict(A=30.5, P=22, C=1.4, g=2.5, dAM=20, Tm=1.2, TM=4, Tw=1.2, TW=1,
             dA=0.16, dP=0.1, dm=1, dM=2, dw=1.3, dW=2.3,
             kA=0.5, kWA=0.14, kMP=0.15, kPm=1, kPw=1, R=2, Smax=0.2, gamma=9)
P2017 = dict(P2025, Tm=1.0, Tw=1.0)

def rhs(t, y, par, S_of_t):
    ATM, P53, mdm2, Mdm2, wip1, Wip1 = np.maximum(y, 0.0)
    S = S_of_t(t)
    dATM = par["A"] * (ATM**2 / (1 + ATM**2 / par["kA"]**2)) * (1 / (1 + Wip1 / par["kWA"])) \
        - par["dA"] * ATM - par["P"] * ATM * Wip1 + S
    dP53 = par["C"] - par["dP"] * P53 - par["g"] * Mdm2 * P53 / (par["kMP"] + P53) / (1 + par["R"] * ATM)
    dmdm2 = par["Tm"] * P53 / (par["kPm"] + P53) - par["dm"] * mdm2
    dMdm2 = par["TM"] * mdm2 - par["dM"] * Mdm2 - par["dAM"] * ATM * Mdm2
    dwip1 = par["Tw"] * P53 / (par["kPw"] + P53) - par["dw"] * wip1
    dWip1 = par["TW"] * wip1 - par["dW"] * Wip1
    return [dATM, dP53, dmdm2, dMdm2, dwip1, dWip1]

def sim(par, S_of_t, T, y0, dt=0.02):
    t_eval = np.arange(0, T + dt, dt)
    sol = solve_ivp(rhs, (0, T), y0, args=(par, S_of_t),
                    method="LSODA", t_eval=t_eval, rtol=1e-7, atol=1e-10, max_step=0.1)
    return sol.t, sol.y

t, y = sim(P2025, lambda t: 0.0, 200.0, [0]*6, dt=0.5)
y0 = y[:, -1]
print("ground state @S=0:", np.round(y0, 4))

print("--- constant S sweep (2025) ---")
for S in [0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.8, 1.2, 2.0, 4.0]:
    t, y = sim(P2025, lambda t, S=S: S, 96.0, y0)
    p53 = y[1][t > 48]
    rel = (np.percentile(p53, 95) - np.percentile(p53, 5)) / (np.median(p53) + 1e-12)
    pk, _ = find_peaks(y[1][t > 48], prominence=0.05)
    ipi = np.diff(t[t > 48][pk]) if len(pk) > 1 else []
    print(f"S={S:>4}: mean={p53.mean():.3f} rel_osc={rel:.3f} n_pk={len(pk)}",
          f"T={np.mean(ipi):.2f}" if len(ipi) else "")

print("--- constant S sweep (2017) ---")
for S in [0.1, 0.2, 0.5, 1.0, 2.0]:
    t, y = sim(P2017, lambda t, S=S: S, 120.0, y0)
    p53 = y[1][t > 48]
    rel = (np.percentile(p53, 95) - np.percentile(p53, 5)) / (np.median(p53) + 1e-12)
    pk, _ = find_peaks(y[1][t > 48], prominence=0.05)
    ipi = np.diff(t[t > 48][pk]) if len(pk) > 1 else []
    print(f"S={S:>4}: mean={p53.mean():.3f} rel_osc={rel:.3f} n_pk={len(pk)}",
          f"T={np.mean(ipi):.2f}" if len(ipi) else "")

print("--- NCS-like decaying DSB input (sat), 2025 ---")
def dsb_ncs(b_s, b_b=2.3, r=0.315, ts=1.0):
    d0 = b_b / r
    d1 = b_s / r + (d0 - b_s / r) * np.exp(-r * ts)
    def f(t):
        if t < ts:
            return b_s / r + (d0 - b_s / r) * np.exp(-r * t)
        return b_b / r + (d1 - b_b / r) * np.exp(-r * (t - ts))
    return f
def Ssat(par, d):
    return par["Smax"] * d / (par["gamma"] + d)
for b_s in [50, 200, 600]:
    f = dsb_ncs(b_s)
    t, y = sim(P2025, lambda t: Ssat(P2025, f(t)), 48.0, y0)
    p53 = y[1]
    p10, p99 = np.percentile(p53, [10, 99])
    pk, props = find_peaks(p53, prominence=0.15 * (p99 - p10), distance=int(2 / 0.02))
    print(f"b_s={b_s}: n_pulse={len(pk)} t_pk={np.round(t[pk],2)} amp={np.round(p53[pk],3)}")
