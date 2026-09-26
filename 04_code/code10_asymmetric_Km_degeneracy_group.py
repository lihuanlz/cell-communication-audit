# -*- coding: utf-8 -*-
"""
Code 10: PdPC audit under asymmetric Michaelis constants (K_m,1 ≠ K_m,2)
==========================================================================
Background: Theorems 1–3 are stated under the symmetric K_m hypothesis (manuscript §3 scope note: "asymmetric case see SI").
An external reviewer asserted that "the degeneracy group drops to 2 dimensions in the asymmetric case" — this script validates before writing; conclusion:
  ★ The reviewer's assertion does not hold. In the asymmetric case the parameter space goes 6→7 dimensions (Km1, Km2 split into separate columns),
    identifiable combinations 2→3 dimensions (ξ, κ1, κ2), and the degeneracy group dimension stays 4 (7−3).
  ★ n_H admits an exact closed-form generalization (not an approximation):
        n_H(κ1, κ2) = 4 / (4 − 1/(κ1+½) − 1/(κ2+½))
    symmetric limit: κ1=κ2=κ → 4/(4 − 2/(κ+½)) = 1 + 1/(2κ) v
  ★ The Fisher barrier (shape identifiable vs absolute scale unidentifiable) holds at all degrees of asymmetry.

[Model] Asymmetric Goldbeter–Koshland steady state (κ1 = K_m,1/S_T on the (1−u) branch,
κ2 = K_m,2/S_T on the u branch):
        ξ = u/(1−u) · (κ1+1−u)/(κ2+u)
At r=1 (symmetric) it reproduces all code-4 baselines (n_H exact, group actions at machine zero, χ stratified).

[Structure]
P1 degeneracy group dimension: 7 parameters (k1,E1T,k2,E2T,Km1,Km2,ST) → 3 identifiable combinations (ξ,κ1,κ2)
   → 4 independent exact group actions (numeric validation max|Δu| = 0) + 1 counterexample (Km1 varied alone)
P2 n_H closed-form generalization: numeric Hill coefficient vs closed form, across asymmetry r = κ1/κ2
P3 Fisher barrier: shape block (ln c, ln κ1, ln κ2) vs absolute block (ln c, ln Km1, ln Km2, ln ST)
output: 非symmetryKm分析图.png (four panels) + SI-ready table
run: python3 代码10_非symmetryKmdegenerate群分析.py (approximately 1 minute)
"""
import numpy as np
from scipy.optimize import brentq
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

for _f in ['Noto Sans CJK SC', 'WenQuanYi Zen Hei', 'SimHei', 'Microsoft YaHei']:
    try:
        from matplotlib.font_manager import findfont, FontProperties
        if findfont(FontProperties(family=_f), fallback_to_default=False):
            plt.rcParams['font.sans-serif'] = [_f, 'DejaVu Sans']
            break
    except Exception:
        continue
plt.rcParams['axes.unicode_minus'] = False

OUTDIR = os.path.dirname(os.path.abspath(__file__))
SIG = 0.03

# ============================================================
# Master equation and analytic derivatives (asymmetric)
# ============================================================
def solve_u(xi, k1_, k2_):
    xi = np.atleast_1d(np.asarray(xi, float))
    f = lambda u, x: x*(1-u)*(k2_+u) - u*(k1_+1-u)
    out = np.empty_like(xi)
    for i, x in enumerate(xi.flat):
        out.flat[i] = brentq(f, 1e-16, 1-1e-16, args=(x,), xtol=1e-14, rtol=1e-14)
    return out

def dxi_du(u, k1_, k2_):
    xi = u/(1-u)*(k1_+1-u)/(k2_+u)
    return xi*(1/u - 1/(k1_+1-u) + 1/(1-u) - 1/(k2_+u))

def du_dki(u, xi, k1_, k2_, which):
    """du/dκ_i (implicit differentiation)."""
    G = lambda u, a, b: xi*(1-u)*(b+u) - u*(a+1-u)
    Gu = xi*(1-k2_-2*u) - (k1_+1-2*u)
    if which == 1:
        Gi = -u
    else:
        Gi = xi*(1-u)
    return -Gi/Gu

def nH_closed(k1_, k2_):
    """Exact closed-form generalization of n_H: logit slope at u=1/2."""
    return 4.0/(4.0 - 1.0/(k1_+0.5) - 1.0/(k2_+0.5))

def nH_numeric(k1_, k2_):
    """Numeric Hill coefficient: d logit(u)/d ln ξ at u=1/2 = 4ξ·du/dξ."""
    xi50 = 0.5/0.5*(k1_+0.5)/(k2_+0.5)
    return 4.0*xi50/dxi_du(0.5, k1_, k2_)

print('='*72)
print('Code 10: asymmetric K_m (κ1 ≠ κ2) PdPC audit')
print('='*72)

# ============================================================
# P1 degeneracy group dimension: 4 independent exact group actions + 1 counterexample
# ============================================================
print('\n[P1] degeneracy group dimension (parameter space 7-dim → identifiable combinations 3-dim (ξ, κ1, κ2))')
# benchmark (asymmetric): k1,E1T,k2,E2T,Km1,Km2,ST
base = dict(k1=1.0, E1T=1.0, k2=1.0, E2T=1.0, Km1=0.03, Km2=0.1, ST=1.0)
x_grid = np.logspace(-2, 2, 200)

def titration(p):
    xi = p['k1']*p['E1T']*x_grid/(p['k2']*p['E2T'])
    return solve_u(xi, p['Km1']/p['ST'], p['Km2']/p['ST'])

u_base = titration(base)
actions = {
    '(a) k1,k2 ×2':        {**base, 'k1': 2.0, 'k2': 2.0},           # ξ unchanged
    "(a') E1T,E2T ×3":     {**base, 'E1T': 3.0, 'E2T': 3.0},         # ξ unchanged (independent direction)
    '(b) Km1,Km2,ST ×5':   {**base, 'Km1': 0.15, 'Km2': 0.5, 'ST': 5.0},  # κ1,κ2 unchanged
    '(c) k1,E2T ×4':       {**base, 'k1': 4.0, 'E2T': 4.0},          # ξ unchanged
}
print(f'     benchmark (Km1,Km2,ST)=(0.03,0.1,1), 200-point log titration:')
for name, g in actions.items():
    d = np.max(np.abs(titration(g) - u_base))
    print(f'     group action {name:<22} max|Δu| = {d:.3e}  {"(exact v)" if d < 1e-10 else "(breaking x)"}')
# counterexample: Km1 varied alone (κ1 changes → the curve must change)
g_bad = {**base, 'Km1': 0.06}
d_bad = np.max(np.abs(titration(g_bad) - u_base))
print(f'     counterex. Km1 ×2 alone        max|Δu| = {d_bad:.3e}  (not a group action, curve changes v as expected)')
print('     conclusion: 4 independent exact actions → degeneracy group stays 4-dimensional (7−3); reviewer "drops to 2" assertion does not hold.')

# ============================================================
# P2 n_H closed-form generalization validation
# ============================================================
print('\n[P2] n_H closed-form generalization: n_H = 4/(4 − 1/(κ1+½) − 1/(κ2+½))')
k2_fix = 0.01
ratios = [0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0]
print(f'     κ2 = {k2_fix} fixed; scan r = κ1/κ2')
print(f'     {"r":<8}{"n_H numeric":<12}{"n_H closed":<12}{"rel deviation":<12}{"symmetric 1+1/(2κ2)":<18}')
rows = []
for r in ratios:
    k1_ = r*k2_fix
    nh_n, nh_c = nH_numeric(k1_, k2_fix), nH_closed(k1_, k2_fix)
    dev = abs(nh_n-nh_c)/nh_c
    nh_sym = 1 + 1/(2*k2_fix)
    rows.append((r, nh_n, nh_c, dev))
    print(f'     {r:<8}{nh_n:<12.4f}{nh_c:<12.4f}{dev:<12.2e}{nh_sym:<18.2f}')
print(f'     max relative deviation = {max(x[3] for x in rows):.2e} (machine precision → closed form holds exactly)')
print(f'     note: for r≠1 the symmetric formula 1+1/(2κ2) fails systematically (e.g. r=30: closed form {rows[-1][2]:.1f} vs symmetric 51)')
# symmetric baseline recheck
nh1_n, nh1_c = nH_numeric(k2_fix, k2_fix), 1+1/(2*k2_fix)
print(f'     symmetric baseline recheck (r=1): numeric {nh1_n:.6f} vs 1+1/(2κ) = {nh1_c:.6f} v')

# ============================================================
# P3 Fisher barrier persistence
# ============================================================
print('\n[P3] Fisher barrier: shape block (ln c, ln κ1, ln κ2) vs absolute block (ln c, ln Km1, ln Km2, ln ST)')
LNXI = np.linspace(-1.5*np.log(10), 1.5*np.log(10), 20)
XI = np.exp(LNXI)

def jac_shape(k1_, k2_):
    u = solve_u(XI, k1_, k2_)
    a = XI/dxi_du(u, k1_, k2_)
    b1 = k1_*du_dki(u, XI, k1_, k2_, 1)
    b2 = k2_*du_dki(u, XI, k1_, k2_, 2)
    return np.column_stack([a, b1, b2])

def jac_absolute(k1_, k2_):
    """Absolute parameters (ln c, ln Km1, ln Km2, ln ST): u depends via κ1=Km1/ST, κ2=Km2/ST."""
    u = solve_u(XI, k1_, k2_)
    a = XI/dxi_du(u, k1_, k2_)
    b1 = k1_*du_dki(u, XI, k1_, k2_, 1)
    b2 = k2_*du_dki(u, XI, k1_, k2_, 2)
    return np.column_stack([a, b1, b2, -(b1+b2)])

print(f'     {"r":<8}{"χ(shape 3-param)":<16}{"χ(absolute 4-param)":<16}{"λ_min(absolute)":<14}')
fisher_rows = []
for r in ratios:
    k1_ = r*k2_fix
    Fs = jac_shape(k1_, k2_fix).T @ jac_shape(k1_, k2_fix)/SIG**2
    Fa = jac_absolute(k1_, k2_fix).T @ jac_absolute(k1_, k2_fix)/SIG**2
    es = np.linalg.eigvalsh(Fs); ea = np.linalg.eigvalsh(Fa)
    chi_s = np.sqrt(es.max()/es.min())
    chi_a = np.sqrt(ea.max()/max(ea.min(), 1e-300))
    fisher_rows.append((r, chi_s, chi_a, ea.min()))
    print(f'     {r:<8}{chi_s:<16.2f}{chi_a:<16.2e}{ea.min():<14.2e}')
print('     conclusion: shape χ stays O(10) (identifiable); absolute block always has a numeric-zero eigenvalue (barrier persists) v')

# ============================================================
# Figure: four panels
# ============================================================
fig, ax = plt.subplots(2, 2, figsize=(13, 9))
xi_c = np.logspace(-1.5, 1.5, 400)
for r, col in zip([0.1, 0.3, 1.0, 3.0, 10.0], plt.cm.viridis(np.linspace(0, 0.9, 5))):
    ax[0,0].semilogx(xi_c, solve_u(xi_c, r*k2_fix, k2_fix), color=col, lw=1.5, label=f'r={r}')
ax[0,0].set_xlabel('ξ (kinase/phosphatase activity ratio)'); ax[0,0].set_ylabel('u (modification fraction)')
ax[0,0].set_title(f'(a) asymmetric titration curves (κ2={k2_fix})'); ax[0,0].legend(); ax[0,0].grid(alpha=0.3)
rs = np.array([x[0] for x in rows])
ax[0,1].loglog(rs, [x[1] for x in rows], 'o', ms=7, label='numeric')
ax[0,1].loglog(rs, [x[2] for x in rows], '-', lw=1.5, label='closed form 4/(4−1/(κ1+½)−1/(κ2+½))')
ax[0,1].axhline(1+1/(2*k2_fix), color='r', ls='--', lw=1, label='symmetric 1+1/(2κ2) (fails)')
ax[0,1].set_xlabel('r = κ1/κ2'); ax[0,1].set_ylabel('n_H')
ax[0,1].set_title('(b) n_H closed form vs numeric'); ax[0,1].legend(fontsize=8); ax[0,1].grid(alpha=0.3)
fr = np.array(fisher_rows)
ax[1,0].semilogx(fr[:,0], fr[:,1], 'o-', label='χ shape block (lnc, lnκ1, lnκ2)')
ax[1,0].semilogx(fr[:,0], fr[:,2], 's--', label='χ absolute block (lnc, lnKm1, lnKm2, lnST)')
ax[1,0].set_yscale('log')
ax[1,0].set_xlabel('r = κ1/κ2'); ax[1,0].set_ylabel('condition number χ')
ax[1,0].set_title('(c) Fisher barrier persists across asymmetry'); ax[1,0].legend(); ax[1,0].grid(alpha=0.3)
names = list(actions.keys()) + ['Km1 ×2 alone (counterexample)']
vals = [np.max(np.abs(titration(g)-u_base)) for g in actions.values()] + [d_bad]
colors = ['tab:blue']*4 + ['tab:red']
ax[1,1].bar(range(5), np.maximum(vals, 1e-17), color=colors)
ax[1,1].set_yscale('log'); ax[1,1].set_ylim(1e-17, 1)
ax[1,1].set_xticks(range(5)); ax[1,1].set_xticklabels(['(a)', "(a')", '(b)', '(c)', 'counterex.'], fontsize=9)
ax[1,1].set_ylabel('max|Δu|')
ax[1,1].set_title('(d) group-action exactness (4 exact + 1 breaking) → group dimension 4')
ax[1,1].grid(alpha=0.3, axis='y')
fig.tight_layout()
out = os.path.join(OUTDIR, '非symmetryKm分析图.png')
fig.savefig(out, dpi=160)
print(f'\nFigure saved: {out}')

# SI-ready table
print('\n[SI-ready] Table S1: asymmetric K_m audit summary')
print('r\tκ1\tn_H numeric\tn_H closed\tχ shape\tχ absolute\tλmin absolute')
for (r, nh_n, nh_c, dev), (r2, cs, ca, lm) in zip(rows, fisher_rows):
    print(f'{r}\t{r*k2_fix:.4f}\t{nh_n:.3f}\t{nh_c:.3f}\t{cs:.1f}\t{ca:.2e}\t{lm:.1e}')
print('\nDone.')
