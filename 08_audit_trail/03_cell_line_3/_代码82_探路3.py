# -*- coding: utf-8 -*-
"""探路3：细扫恒定 S + 长窗口 + 高初值，验证振荡区是否存在。"""
import numpy as np
from scipy.integrate import solve_ivp

P2025 = dict(A=30.5, P=22, C=1.4, g=2.5, dAM=20, Tm=1.2, TM=4, Tw=1.2, TW=1,
             dA=0.16, dP=0.1, dm=1, dM=2, dw=1.3, dW=2.3,
             kA=0.5, kWA=0.14, kMP=0.15, kPm=1, kPw=1, R=2, Smax=0.2, gamma=9)

def rhs(t, y, par, S):
    ATM, P53, mdm2, Mdm2, wip1, Wip1 = np.maximum(y, 0.0)
    dATM = par["A"] * (ATM**2 / (1 + ATM**2 / par["kA"]**2)) * (1 / (1 + Wip1 / par["kWA"])) \
        - par["dA"] * ATM - par["P"] * ATM * Wip1 + S
    dP53 = par["C"] - par["dP"] * P53 - par["g"] * Mdm2 * P53 / (par["kMP"] + P53) / (1 + par["R"] * ATM)
    dmdm2 = par["Tm"] * P53 / (par["kPm"] + P53) - par["dm"] * mdm2
    dMdm2 = par["TM"] * mdm2 - par["dM"] * Mdm2 - par["dAM"] * ATM * Mdm2
    dwip1 = par["Tw"] * P53 / (par["kPw"] + P53) - par["dw"] * wip1
    dWip1 = par["TW"] * wip1 - par["dW"] * Wip1
    return [dATM, dP53, dmdm2, dMdm2, dwip1, dWip1]

def sim(par, S, T, y0, dt=0.05):
    t_eval = np.arange(0, T + dt, dt)
    sol = solve_ivp(rhs, (0, T), y0, args=(par, S),
                    method="LSODA", t_eval=t_eval, rtol=1e-8, atol=1e-11, max_step=0.05)
    return sol.t, sol.y

y0_low = [0, 0.4363, 0.3645, 0.7291, 0.2804, 0.1219]
y0_high = [3.0, 8.0, 2.0, 1.0, 2.0, 1.0]

for S in np.round(np.arange(0.02, 0.62, 0.04), 3):
    for tag, y0 in [("low", y0_low), ("high", y0_high)]:
        t, y = sim(P2025, S, 300.0, y0, dt=0.25)
        p53 = y[1][t > 200]
        rel = (p53.max() - p53.min()) / (np.median(p53) + 1e-12)
        print(f"S={S:.2f} ic={tag}: mean={p53.mean():.3f} rel_osc={rel:.3f}")
