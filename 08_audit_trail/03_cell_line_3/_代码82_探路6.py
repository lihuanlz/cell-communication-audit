# -*- coding: utf-8 -*-
"""探路6：勘误后的方程（kA 非 kA^2；p53 保护因子 (1+R/(1+ATM*))）验证振荡行为。"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

P2025 = dict(A=30.5, P=22, C=1.4, g=2.5, dAM=20, Tm=1.2, TM=4, Tw=1.2, TW=1,
             dA=0.16, dP=0.1, dm=1, dM=2, dw=1.3, dW=2.3,
             kA=0.5, kWA=0.14, kMP=0.15, kPm=1, kPw=1, R=2, Smax=0.2, gamma=9)
P2017 = dict(P2025, Tm=1.0, Tw=1.0)

def rhs(t, y, par, S):
    ATM, P53, mdm2, Mdm2, wip1, Wip1 = np.maximum(y, 0.0)
    hill = ATM**2 / (1.0 + ATM**2 / par["kA"])          # 勘误：/kA 非 /kA^2
    dATM = par["A"] * hill / (1.0 + Wip1 / par["kWA"]) \
        - par["dA"] * ATM - par["P"] * ATM * Wip1 + S
    dP53 = par["C"] - par["dP"] * P53 \
        - par["g"] * Mdm2 * P53 / (par["kMP"] + P53) * (1.0 + par["R"] / (1.0 + ATM))  # 勘误
    dmdm2 = par["Tm"] * P53 / (par["kPm"] + P53) - par["dm"] * mdm2
    dMdm2 = par["TM"] * mdm2 - par["dM"] * Mdm2 - par["dAM"] * ATM * Mdm2
    dwip1 = par["Tw"] * P53 / (par["kPw"] + P53) - par["dw"] * wip1
    dWip1 = par["TW"] * wip1 - par["dW"] * Wip1
    return [dATM, dP53, dmdm2, dMdm2, dwip1, dWip1]

def sim(par, S, T, y0, dt=0.05):
    t_eval = np.arange(0, T + dt, dt)
    sol = solve_ivp(rhs, (0, T), y0, args=(par, S), method="LSODA",
                    t_eval=t_eval, rtol=1e-8, atol=1e-11, max_step=0.05)
    return sol.t, sol.y

t, y = sim(P2025, 0.0, 300.0, [0]*6, dt=1.0)
y0 = y[:, -1]
print("ground state @S=0:", np.round(y0, 4))

print("--- constant S sweep (2025 params) ---")
for S in [0.02, 0.05, 0.08, 0.1, 0.15, 0.2, 0.3, 0.5, 1.0]:
    t, y = sim(P2025, S, 200.0, y0)
    p53 = y[1][t > 100]
    rel = (p53.max() - p53.min()) / (np.median(p53) + 1e-12)
    pk, _ = find_peaks(p53, prominence=0.05)
    ipi = np.diff(t[t > 100][pk]) if len(pk) > 1 else []
    print(f"S={S:5.2f}: mean={p53.mean():.3f} rel_osc={rel:.3f} n_pk={len(pk)}",
          f"T={np.mean(ipi):.2f}+-{np.std(ipi):.2f}" if len(ipi) else "")

print("--- 2017 params ---")
for S in [0.05, 0.1, 0.2]:
    t, y = sim(P2017, S, 200.0, y0)
    p53 = y[1][t > 100]
    pk, _ = find_peaks(p53, prominence=0.05)
    ipi = np.diff(t[t > 100][pk]) if len(pk) > 1 else []
    print(f"S={S:5.2f}: n_pk={len(pk)}",
          f"T={np.mean(ipi):.2f}" if len(ipi) else "")
