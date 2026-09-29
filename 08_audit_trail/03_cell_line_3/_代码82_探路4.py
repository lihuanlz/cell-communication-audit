# -*- coding: utf-8 -*-
"""探路4：不动点 + Jacobian 特征值随 S 的稳定性扫描，定位 Hopf/SNIC。"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

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
    return np.array([dATM, dP53, dmdm2, dMdm2, dwip1, dWip1])

def fp(par, S, y0):
    sol = solve_ivp(rhs, (0, 2000), y0, args=(par, S), method="LSODA",
                    rtol=1e-10, atol=1e-13)
    return sol.y[:, -1]

def jac(par, S, y, h=1e-6):
    n = len(y)
    J = np.zeros((n, n))
    f0 = rhs(0, y, par, S)
    for i in range(n):
        dy = np.zeros(n); dy[i] = h * max(abs(y[i]), 1e-3)
        J[:, i] = (rhs(0, y + dy, par, S) - rhs(0, y - dy, par, S)) / (2 * dy[i])
    return J

y0 = [0, 0.4363, 0.3645, 0.7291, 0.2804, 0.1219]
print(" S      maxRe(eig)  imag?  ATM*   p53")
for S in [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 1.0]:
    y = fp(P2025, S, y0)
    ev = np.linalg.eigvals(jac(P2025, S, y))
    i = np.argmax(ev.real)
    print(f"{S:5.2f}  {ev[i].real:+.5f}  {abs(ev[i].imag):.4f}  {y[0]:.4f} {y[1]:.4f}")
