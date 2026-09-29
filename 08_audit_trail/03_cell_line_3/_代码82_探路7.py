# -*- coding: utf-8 -*-
"""探路7：勘误版完整行为谱。S 上限扫描（R4 转变点）、NCS 脉冲串、随机基态。"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

P2025 = dict(A=30.5, P=22, C=1.4, g=2.5, dAM=20, Tm=1.2, TM=4, Tw=1.2, TW=1,
             dA=0.16, dP=0.1, dm=1, dM=2, dw=1.3, dW=2.3,
             kA=0.5, kWA=0.14, kMP=0.15, kPm=1, kPw=1, R=2, Smax=0.2, gamma=9)

def rhs(t, y, par, S_of_t):
    ATM, P53, mdm2, Mdm2, wip1, Wip1 = np.maximum(y, 0.0)
    S = S_of_t(t)
    hill = ATM**2 / (1.0 + ATM**2 / par["kA"])
    dATM = par["A"] * hill / (1.0 + Wip1 / par["kWA"]) \
        - par["dA"] * ATM - par["P"] * ATM * Wip1 + S
    dP53 = par["C"] - par["dP"] * P53 \
        - par["g"] * Mdm2 * P53 / (par["kMP"] + P53) * (1.0 + par["R"] / (1.0 + ATM))
    dmdm2 = par["Tm"] * P53 / (par["kPm"] + P53) - par["dm"] * mdm2
    dMdm2 = par["TM"] * mdm2 - par["dM"] * Mdm2 - par["dAM"] * ATM * Mdm2
    dwip1 = par["Tw"] * P53 / (par["kPw"] + P53) - par["dw"] * wip1
    dWip1 = par["TW"] * wip1 - par["dW"] * Wip1
    return [dATM, dP53, dmdm2, dMdm2, dwip1, dWip1]

def sim(par, S_of_t, T, y0, dt=0.05, max_step=0.05):
    t_eval = np.arange(0, T + dt, dt)
    sol = solve_ivp(rhs, (0, T), y0, args=(par, S_of_t), method="LSODA",
                    t_eval=t_eval, rtol=1e-7, atol=1e-10, max_step=max_step)
    return sol.t, sol.y

t, y = sim(P2025, lambda t: 0.0, 300.0, [0]*6, dt=1.0)
y0 = y[:, -1]

print("--- high-S: oscillation death (R4 transition point) ---")
for S in [0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0]:
    t, y = sim(P2025, lambda t, S=S: S, 200.0, y0)
    p53 = y[1][t > 100]
    rel = (p53.max() - p53.min()) / (np.median(p53) + 1e-12)
    pk, _ = find_peaks(p53, prominence=0.05)
    ipi = np.diff(t[t > 100][pk]) if len(pk) > 1 else []
    print(f"S={S:4.1f}: mean={p53.mean():.3f} rel_osc={rel:.3f} n_pk={len(pk)}",
          f"T={np.mean(ipi):.2f}" if len(ipi) else "")

print("--- NCS deterministic (sat), b_s sweep: pulse trains ---")
def dsb_ncs(b_s, b_b=2.3, r=0.315, ts=1.0):
    d0 = b_b / r
    d1 = b_s / r + (d0 - b_s / r) * np.exp(-r * ts)
    def f(t):
        if t < ts:
            return b_s / r + (d0 - b_s / r) * np.exp(-r * t)
        return b_b / r + (d1 - b_b / r) * np.exp(-r * (t - ts))
    return f
def Ssat(d):
    return P2025["Smax"] * d / (P2025["gamma"] + d)
for b_s in [50, 200, 600]:
    f = dsb_ncs(b_s)
    t, y = sim(P2025, lambda t: Ssat(f(t)), 48.0, y0)
    p53 = y[1]
    p10, p99 = np.percentile(p53, [10, 99])
    pk, _ = find_peaks(p53, prominence=0.15 * (p99 - p10), distance=int(2 / 0.05))
    print(f"b_s={b_s}: n_pulse={len(pk)} t={np.round(t[pk],1)} amp={np.round(p53[pk],2)}")

print("--- NCS deterministic (log), b_s sweep ---")
def Slog(d):
    return P2025["Smax"] * np.log(d / P2025["gamma"] + 1.0)
for b_s in [50, 200, 600]:
    f = dsb_ncs(b_s)
    t, y = sim(P2025, lambda t: Slog(f(t)), 48.0, y0)
    p53 = y[1]
    p10, p99 = np.percentile(p53, [10, 99])
    pk, _ = find_peaks(p53, prominence=0.15 * (p99 - p10), distance=int(2 / 0.05))
    print(f"b_s={b_s}: n_pulse={len(pk)} t={np.round(t[pk],1)} amp={np.round(p53[pk],2)} "
          f"late_mean={p53[t>36].mean():.2f}")

print("--- stochastic basal DSB (two-pool Gillespie), 5 cells ---")
RC = np.log(2) / 20
def dsb_stoch_grid(rng, t_grid, b_s=0.0, b_basal=0.7, t_switch=1.0, r_i=0.35, r_c=RC, f_c=0.1):
    dt = t_grid[1] - t_grid[0]
    out = np.zeros_like(t_grid)
    n_i = rng.poisson(b_basal * (1 - f_c) / r_i)
    n_c = rng.poisson(b_basal * f_c / r_c)
    t = 0.0; gi = 0
    T_end = t_grid[-1]
    while t < T_end:
        b = b_s if t < t_switch else b_basal
        a0 = b + r_i * n_i + r_c * n_c
        tau = rng.exponential(1.0 / a0) if a0 > 0 else np.inf
        t_next = t + tau
        t_step_end = min(t_next, t_switch if t < t_switch else T_end, T_end)
        while gi < len(t_grid) and t_grid[gi] < t_step_end:
            out[gi] = n_i + n_c; gi += 1
        if t_next > t_step_end - 1e-12 and t_step_end < T_end:
            t = t_step_end
            continue
        if t_next >= T_end:
            break
        u = rng.random() * a0
        if u < b:
            if rng.random() < f_c: n_c += 1
            else: n_i += 1
        elif u < b + r_i * n_i:
            n_i -= 1
        else:
            n_c -= 1
        t = t_next
    while gi < len(t_grid):
        out[gi] = n_i + n_c; gi += 1
    return out

rng = np.random.default_rng(0)
t_grid = np.arange(0, 48.0001, 0.05)
for cell in range(5):
    grid = dsb_stoch_grid(rng, t_grid)
    def f(t, g=grid, tg=t_grid):
        i = int(t / 0.05)
        return g[min(i, len(g) - 1)]
    t, y = sim(P2025, lambda tt: Ssat(f(tt)), 48.0, y0)
    p53 = y[1]
    pk, _ = find_peaks(p53, prominence=0.3, height=0.8, distance=int(2 / 0.05))
    print(f"cell{cell}: maxDSB={grid.max():.0f} n_pulse={len(pk)} t={np.round(t[pk],1)} amp={np.round(p53[pk],2)}")
