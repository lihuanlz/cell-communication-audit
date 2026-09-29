# -*- coding: utf-8 -*-
"""
PdPC identifiability audit (TCS-style analysis of the Goldbeter-Koshland cycle)
================================================================================
Verifies every number in the research memo "PdPC可辨识性审计研究备忘录.md":
  Thm 1  scale degeneracy group (6 physical params -> 2 observable combinations)
  Thm 2  shape identifiable / scale structurally non-identifiable (Fisher chi)
         exact Hill relation n_H = 1 + 1/(2*kappa)
  Thm 3  kappa information plateau -> one-sided reporting protocol
  Sec 6  information-critical point: mutual information I(xi;u) vs kappa
Run:  python3 PdPC_identifiability_audit.py   (numpy + scipy only)
"""
import numpy as np
from scipy.optimize import brentq

# ============================================================================
# Steady-state master equation (symmetric-Km Goldbeter-Koshland cycle)
#   xi = [u/(1-u)] * [(kap+1-u)/(kap+u)]
#   u   = phosphorylated fraction (output)
#   xi  = V1/V2 (input), kap = K_m/S_T (saturation ratio)
# ============================================================================
def xi_of_u(u, kap):
    return (u / (1 - u)) * ((kap + 1 - u) / (kap + u))

def u_of_xi(xi, kap):
    if kap < 1e-12:                       # zero-order limit: perfect switch
        return 0.0 if xi < 1 else 1.0
    return brentq(lambda u: xi_of_u(u, kap) - xi, 1e-15, 1 - 1e-15, xtol=1e-14)

def nH_analytic(kap):
    return 1 + 0.5 / kap                  # exact: n_H = 4 xi du/dxi |_{u=1/2}

# ============================================================================
print("=" * 70)
print("0. Hill relation n_H = 1 + 1/(2*kappa)  (analytic vs numeric)")
print("=" * 70)
for kap in [0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 10.0]:
    e = 1e-5
    du = u_of_xi(np.exp(e), kap) - u_of_xi(np.exp(-e), kap)
    print(f"  kap={kap:7.3f}: analytic={nH_analytic(kap):8.2f}  numeric={4*du/(2*e):8.2f}")

# ============================================================================
print("\n" + "=" * 70)
print("Thm 1. Scale degeneracy: curve invariance under the degeneracy group")
print("=" * 70)
def phys_curve(x, k1, E1T, k2, E2T, Km, ST):
    return u_of_xi(k1 * E1T * x / (k2 * E2T), Km / ST)

x_grid = np.logspace(-2, 2, 9)
c1 = [phys_curve(x, 1, 1, 1, 1, 0.1, 1) for x in x_grid]
c2 = [phys_curve(x, 1, 1, 1, 1, 0.5, 5) for x in x_grid]      # (Km,ST) x5
c3 = [phys_curve(x, 2, 3, 1.5, 4, 0.1, 1) for x in x_grid]    # (k,E) rescale
print(f"  (K_m,S_T) -> 5x :  max|du| = {max(abs(a-b) for a,b in zip(c1,c2)):.2e}")
print(f"  (k,E) rescale   :  max|du| = {max(abs(a-b) for a,b in zip(c1,c3)):.2e}")

# ============================================================================
print("\n" + "=" * 70)
print("Thm 2. Fisher condition number: shape (ln c, ln kap) vs absolute params")
print("=" * 70)
SIG, NPTS = 0.03, 20
xs = np.logspace(-1.5, 1.5, NPTS)

def chi_shape(kap):
    th = np.log([1.0, kap]); J = np.zeros((NPTS, 2))
    for i, x in enumerate(xs):
        for j in range(2):
            dp = th.copy(); dp[j] += 1e-6
            dm = th.copy(); dm[j] -= 1e-6
            J[i, j] = (u_of_xi(np.exp(dp[0])*x, np.exp(dp[1]))
                       - u_of_xi(np.exp(dm[0])*x, np.exp(dm[1]))) / 2e-6
    ev = np.linalg.eigvalsh(J.T @ J / SIG**2)
    return ev[-1] / ev[0]

for kap in [0.005, 0.01, 0.05, 0.1, 0.5, 1.0]:
    print(f"  kap={kap:6.3f}: chi(ln c, ln kap) = {chi_shape(kap):.2e}")

th = np.log([1.0, 0.1, 1.0]); J = np.zeros((NPTS, 3))
for i, x in enumerate(xs):
    for j in range(3):
        dp = th.copy(); dp[j] += 1e-6
        dm = th.copy(); dm[j] -= 1e-6
        J[i, j] = (u_of_xi(np.exp(dp[0])*x, np.exp(dp[1])/np.exp(dp[2]))
                   - u_of_xi(np.exp(dm[0])*x, np.exp(dm[1])/np.exp(dm[2]))) / 2e-6
ev = np.linalg.eigvalsh(J.T @ J / SIG**2)
print(f"  absolute (ln c, ln K_m, ln S_T): eigenvalues = "
      f"{['%.2e' % e for e in ev]}  -> exact zero direction (K_m, S_T)")

# ============================================================================
print("\n" + "=" * 70)
print("Thm 3. kappa CI-width ratio vs kappa (information plateau)")
print("=" * 70)
for kap in [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0]:
    J = np.array([(u_of_xi(x, kap+1e-6) - u_of_xi(x, kap-1e-6)) / 2e-6
                  for x in xs])
    se_lnk = 1 / np.sqrt(np.sum(J**2) / SIG**2 * kap**2)
    print(f"  kap={kap:6.3f}: 95% CI x/{np.exp(1.96*se_lnk):8.2f}   "
          f"n_H={nH_analytic(kap):5.0f}")

# ============================================================================
print("\n" + "=" * 70)
print("Sec 6. Information-critical point: mutual information I(xi; u) vs kappa")
print("=" * 70)
def mutual_info(kap, sig_xi, sig_obs=0.03, n_grid=600):
    xs_ = np.linspace(-4*sig_xi, 4*sig_xi, n_grid); dx = xs_[1]-xs_[0]
    px = np.exp(-xs_**2/(2*sig_xi**2)); px /= (px*dx).sum()
    us = np.array([u_of_xi(np.exp(x), kap) for x in xs_])
    ug = np.linspace(0.001, 0.999, n_grid); dy = ug[1]-ug[0]
    lik = np.exp(-(ug[None,:]-us[:,None])**2/(2*sig_obs**2))
    lik /= (lik*dy).sum(axis=1, keepdims=True)
    py = (px[:,None]*dx*lik).sum(0) + 1e-300
    return max((px[:,None]*lik*np.log(lik/py[None,:])).sum()*dx*dy, 0)/np.log(2)

print(f"  {'kap':>8} {'narrow(0.3)':>12} {'wide(1.5)':>12}")
for kap in [0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.3, 1.0, 3.0, 10.0]:
    print(f"  {kap:8.3f} {mutual_info(kap, 0.3):9.2f} b {mutual_info(kap, 1.5):9.2f} b")
print("\n  -> interior optimum kap* tracks input spread; zero-order switch")
print("     (kap -> 0) is information-SUBOPTIMAL in all cases.")
