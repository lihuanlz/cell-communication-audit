# -*- coding: utf-8 -*-
"""
Direction 2: Cooperative critical phenomena in the canonical TCS ensemble
=========================================================================
Verifies the rigorous critical-cooperativity analysis:

  (1) NO-GO THEOREM (any finite binding polynomial Phi, e.g. MWC):
        xi'(p) = n x / Var(j) + 1/kappa > 0   always
      -> no binding phase transition at any kappa (variance identity).

  (2) Bragg-Williams mean-field attraction (infinite-range cooperative unit):
        generalized master equation  xi(p) = [p/(1-p)] exp(-2cp) + p/kappa
        spinodal:  z'(p) = -1/kappa,   z'(p) = exp(-2cp)[1 - 2cp(1-p)]/(1-p)^2
        dip depth g(c) = -min_p z'(p;c) is bounded:  g(c) -> e^-2  as c -> infinity
      => ABSOLUTE FLOOR: no transition at any c for kappa <= e^2 = 7.389
      => for kappa > e^2: c*(kappa) solves g(c*) = 1/kappa
         asymptotics: c* = 2 + e^2/(2 kappa)      (kappa -> infinity)
                      c* ~ sqrt(kappa / 2(kappa - e^2))  (kappa -> e^2+)
         p_c: 1/2 -> 2/c* -> 0.

  (3) Canonical first-order transition at (kappa=10, c=4): three equilibria,
      spinodals, Maxwell tie point m*, occupancy jump.

  (4) Depletion clamps the apparent Hill divergence: n_H^max(kappa) finite
      below kappa = 1/g(c), diverges exactly at the spinodal threshold.

All equations follow the document's notation:
  xi = L_T/K (dimensionless total feed), kappa = K/(n R_T),
  generalized TCS master equation (S10.21): f xi = x + nu(x)/(n kappa),
  Landau free energy (S1d.6): F = -ln(1-p) - p + p^2/(2 kappa), xi = dF/dp.
"""

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.integrate import quad

E2 = np.e ** 2          # kappa_c
E_2 = np.e ** (-2)      # dip-depth bound


# ----------------------------------------------------------------------
# Bragg-Williams cooperative activity and its derivative
# ----------------------------------------------------------------------
def z_bw(p, c):
    """Cooperative activity z(p) = [p/(1-p)] exp(-2cp)."""
    return (p / (1 - p)) * np.exp(-2 * c * p)


def zprime_bw(p, c):
    """dz/dp = exp(-2cp) [1 - 2cp(1-p)] / (1-p)^2   (checked against FD)."""
    return np.exp(-2 * c * p) * (1 - 2 * c * p * (1 - p)) / (1 - p) ** 2


def xi(p, c, kap):
    """Canonical TCS isotherm (generalized master equation)."""
    return z_bw(p, c) + p / kap


def xi1(p, c, kap):
    return zprime_bw(p, c) + 1.0 / kap


# ----------------------------------------------------------------------
# (2) dip depth g(c) and critical cooperativity c*(kappa)
# ----------------------------------------------------------------------
def gmin(c):
    """g(c) = -min_p z'(p;c), and the argmin. Negative-dip interval
    2cp(1-p) > 1 is bracketed analytically for robustness."""
    if c <= 2:
        return 0.0, 0.5
    disc = np.sqrt(1 - 2.0 / c)
    lo, hi = (1 - disc) / 2, (1 + disc) / 2
    grid = np.linspace(lo, hi, 4000)
    i = int(np.argmin(zprime_bw(grid, c)))
    a, b = grid[max(i - 2, 0)], grid[min(i + 2, len(grid) - 1)]
    res = minimize_scalar(lambda p: zprime_bw(p, c), bounds=(a, b),
                          method='bounded', options={'xatol': 1e-16})
    return -res.fun, res.x


def c_star(kap):
    """Critical cooperativity: g(c*) = 1/kappa, defined for kappa > e^2."""
    return brentq(lambda c: gmin(c)[0] - 1.0 / kap, 2.0 + 1e-9, 1e5,
                  xtol=1e-11, rtol=1e-12)


# ----------------------------------------------------------------------
# (3) canonical equilibria and Maxwell construction
# ----------------------------------------------------------------------
def equilibria(m, kap, c):
    """Solve z_bw(p) = (m - p)/kappa  (equilibrium: xi EOS + depletion identity)."""
    f = lambda p: z_bw(p, c) - (m - p) / kap
    grid = np.linspace(1e-10, min(m, 1 - 1e-10), 20000)
    vals = f(grid)
    roots = []
    for i in range(len(grid) - 1):
        if vals[i] * vals[i + 1] < 0:
            roots.append(brentq(f, grid[i], grid[i + 1], xtol=1e-14))
    return roots


def Fprime(p, m, kap, c):
    """Canonical free-energy derivative F'(p) at fixed feed m (Stirling
    saddle exponent of the TCS partition function + BW attraction)."""
    return np.log(p * np.exp(-2 * c * p) / ((m - p) * (1 - p))) + np.log(kap)


def Fcan(p, m, kap, c):
    val, _ = quad(lambda pp: Fprime(pp, m, kap, c), 1e-10, p, epsabs=1e-10)
    return val


# ----------------------------------------------------------------------
# (4) apparent Hill coefficient
# ----------------------------------------------------------------------
def nH_bw(p, c, kap):
    """n_H = d ln(p/(1-p)) / d ln xi = xi / [p(1-p) xi'(p)]."""
    return xi(p, c, kap) / (p * (1 - p) * xi1(p, c, kap))


def nH_max_bw(c, kap):
    res = minimize_scalar(lambda p: -nH_bw(p, c, kap), bounds=(1e-9, 1 - 1e-9),
                          method='bounded')
    return -res.fun, res.x


# ----------------------------------------------------------------------
# (1) MWC binding polynomial: variance-identity no-go
# ----------------------------------------------------------------------
def mwc_nu(x, n, L, c0):
    """nu(x) for Phi = (1+x)^n + L(1+c0 x)^n."""
    R = (1 + x) ** n
    T = L * (1 + c0 * x) ** n
    return (n * x * (1 + x) ** (n - 1) + L * n * c0 * x * (1 + c0 * x) ** (n - 1)) / (R + T)


if __name__ == "__main__":
    print("=" * 72)
    print("(0) analytic z'(p) vs finite differences")
    ps = np.linspace(0.01, 0.99, 9)
    h = 1e-7
    for c in [1.5, 3.0]:
        fd = (z_bw(ps + h, c) - z_bw(ps - h, c)) / (2 * h)
        an = zprime_bw(ps, c)
        print(f"  c={c}: max rel err = {np.max(np.abs(fd - an) / np.abs(an)):.2e}")

    print("=" * 72)
    print("(2a) dip depth g(c): monotone 0 -> e^-2 (bounded!)")
    for c in [2.01, 2.5, 3, 4, 5, 10, 100, 1e5]:
        g, pm = gmin(c)
        print(f"  c={c:9.2f}: g={g:.8f}  p_min={pm:.6f}  e^-2-g={E_2 - g:.2e}")
    print(f"  e^-2 = {E_2:.6f}  =>  kappa_c = e^2 = {E2:.4f}")

    print("=" * 72)
    print("(2b) c*(kappa) curve with criticality checks xi'=xi''=0, xi'''!=0")
    print(f"{'kappa':>10} {'c*(k)':>10} {'p_c':>9} {'xi1':>9} {'xi3':>10} "
          f"{'2+e^2/2k':>10} {'sqrt-law':>9}")
    for kap in [7.392, 7.5, 8, 10, 20, 50, 100, 1e3, 1e6]:
        cs = c_star(kap)
        g, pc = gmin(cs)
        h2 = 1e-4
        x3 = (zprime_bw(pc + 2 * h2, cs) - 2 * zprime_bw(pc + h2, cs)
              + 2 * zprime_bw(pc - h2, cs) - zprime_bw(pc - 2 * h2, cs)) / (2 * h2 ** 3)
        print(f"{kap:10.3f} {cs:10.4f} {pc:9.5f} {xi1(pc, cs, kap):9.1e} {x3:10.2e} "
              f"{2 + E2 / (2 * kap):10.5f} {np.sqrt(kap / (2 * (kap - E2))):9.4f}")

    print("=" * 72)
    print("(3) canonical first-order transition at kappa=10, c=4")
    kap, c = 10.0, 4.0
    _, pmin4 = gmin(c)
    lo, hi = (1 - np.sqrt(1 - 2 / c)) / 2, (1 + np.sqrt(1 - 2 / c)) / 2
    p_sp = [brentq(lambda p: zprime_bw(p, c) + 1 / kap, lo + 1e-9, pmin4),
            brentq(lambda p: zprime_bw(p, c) + 1 / kap, pmin4, hi - 1e-9)]
    m_sp = [p + kap * z_bw(p, c) for p in p_sp]
    print(f"  spinodal p = {p_sp[0]:.4f}, {p_sp[1]:.4f};  "
          f"bistable m-window = [{min(m_sp):.4f}, {max(m_sp):.4f}]")

    def min_gap(m):
        r = equilibria(m, kap, c)
        if len(r) < 3:
            return np.nan
        return Fcan(r[2], m, kap, c) - Fcan(r[0], m, kap, c)

    m_tie = brentq(min_gap, min(m_sp) + 1e-4, max(m_sp) - 1e-4, xtol=1e-8)
    r = equilibria(m_tie, kap, c)
    print(f"  Maxwell tie m* = {m_tie:.5f};  jump p: {r[0]:.4f} -> {r[2]:.4f}"
          f"  (Delta p = {r[2] - r[0]:.4f})")

    print("=" * 72)
    print("(4) depletion clamps the Hill divergence (c=2.5, supercritical GC)")
    kap_div = 1.0 / gmin(2.5)[0]
    print(f"  threshold kappa = 1/g(2.5) = {kap_div:.3f}")
    for k in [0.5, 1, 2, 5, 8, 10, 11, 11.9, 11.99, 11.999]:
        nHm, pm = nH_max_bw(2.5, k)
        print(f"  kappa={k:7.3f}: n_H^max = {nHm:12.2f} at p={pm:.3f}")
    print("  kappa >= threshold: n_H -> infinity (spinodal; first-order transition)")

    print("=" * 72)
    print("(1) MWC no-go: Var(j) > 0  =>  z'(p) > 0  =>  no spinodal at any kappa")
    xs = np.logspace(-6, 4, 4000)
    n_m, L_m, c0_m = 4, 1e3, 0.05
    hh = 1e-6
    varj = (mwc_nu(xs * np.exp(hh), n_m, L_m, c0_m)
            - mwc_nu(xs * np.exp(-hh), n_m, L_m, c0_m)) / (2 * hh)
    dzdp = n_m * xs / varj
    print(f"  MWC n={n_m}, L={L_m}, c0={c0_m}: min dz/dp = {dzdp.min():.4f} > 0")
    p_arr = mwc_nu(xs, n_m, L_m, c0_m) / n_m
    nH_gc = np.gradient(np.log(p_arr / (1 - p_arr)), np.log(xs))
    print(f"  grand-canonical n_H^max = {nH_gc.max():.3f} (finite; never diverges)")
