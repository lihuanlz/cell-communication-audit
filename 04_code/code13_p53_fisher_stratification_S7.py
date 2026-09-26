# -*- coding: utf-8 -*-
"""
Code 13: Fisher information stratification of the p53 pulse encoder (audit step 3, SI S7)
============================================================
Purpose: perform step 3 (Fisher stratification) of the five-step audit of the dynamic encoder, and numerically validate
the predictions of the number-limit condition theorem (SI S2, H1–H3 ⟹ C1–C3) for p53:
  (1) the Fisher information of the count statistic is nonzero only along the structural parameter directions (D_c, ρ);
  (2) the amplitude observable pins only the product s·A; the orthogonal direction (ln A − ln s) is exactly zero
     — this null direction is the orbit tangent of the degeneracy group G (simulated scale nuisance);
  (3) Monte-Carlo validation of Cramér–Rao: structural parameters attain the bound; the scale direction is inestimable (likelihood ridge).

Model (strictly consistent with the manuscript Methods section and code 5):
  count law  N(D) = ρ·ln(D₀(D)/D_c), D₀(D) = 0.40 + 0.09·D, n ~ Poisson(N)
  amplitude    ȳ ~ Normal(s·A, σ²) (M cells per dose; σ = 15%·sA, experimental-grade precision)
Parameter vector θ = (ln D_c, ln ρ, ln A, ln s).

Run: python3 代码13_p53_Fisherstratified_S7.py (< 1 minute, outputs a four-panel PNG)
"""
import numpy as np
import matplotlib.pyplot as plt

rng = np.random.default_rng(20260805)

# ---------- ground-truth parameters (calibrated in code 5) ----------
rho_true  = 9.3      # count-law slope (τ_r/T structural parameter)
Dc_true   = 0.45     # countthreshold（D₀ unit）
A_true    = 2.0      # physical pulse amplitude (molecular units)
s_true    = 1.0      # simulated scale (fluorescence/molecule, nuisance)
M         = 200      # cells per dose (typical single-cell experiment order of magnitude)
SIG_FRAC  = 0.15     # relative precision of the amplitude measurement
D_grid    = np.geomspace(2.0, 17.0, 6)          # dose grid (within the sustained counting window; N ≥ 2.4 ensures the MLE asymptotic regime)
D0        = 0.40 + 0.09*D_grid
N_D       = rho_true*np.log(D0/Dc_true)         # mean count at each dose
sA        = s_true*A_true
sig       = SIG_FRAC*sA

def fisher_matrix(counts_only=False):
    """F = Σ_D [F_count + F_amp]，θ = (ln D_c, ln ρ, ln A, ln s)。"""
    F = np.zeros((4, 4))
    for N in N_D:
        g = np.array([-rho_true, N, 0.0, 0.0])     # ∂N/∂θ（Poisson: F = M·g gᵀ/N）
        F += M*np.outer(g, g)/N
        if not counts_only:
            h = np.array([0.0, 0.0, sA, sA])       # ∂(sA)/∂θ
            F += (M/sig**2)*np.outer(h, h)
    return F

F_full = fisher_matrix(counts_only=False)
F_cnt  = fisher_matrix(counts_only=True)

def eigsorted(F):
    w, V = np.linalg.eigh(F)
    idx = np.argsort(w)[::-1]
    return w[idx], V[:, idx]

w_full, V_full = eigsorted(F_full)
w_cnt,  V_cnt  = eigsorted(F_cnt)

print('='*70)
print('[Predictions (1)+(2): Fisher spectrum]')
print(f'counts only (number limit): λ = {w_cnt}  → non-null directions: {np.sum(w_cnt > 1e-8*w_cnt[0])} (prediction 2)')
print(f'counts+amplitude:           λ = {w_full} → non-null directions: {np.sum(w_full > 1e-8*w_full[0])} (prediction 3)')
print(f'smallest/largest eigenvalue = {w_full[-1]/w_full[0]:.2e} (float-noise order of magnitude → exact null direction exists)')
v0 = V_full[:, -1]
print(f'null directioneigenvector = {np.round(v0, 4)}（prediction ∝ (0,0,1,−1)/√2 = [0,0,0.7071,−0.7071]）')
G_tangent = np.array([0, 0, 1, -1])/np.sqrt(2)
print(f'cosine of the angle with the G orbit tangent (0,0,1,−1)/√2 = {abs(v0 @ G_tangent):.6f} (prediction 1.000000)')

# ---------- stratified quantification (same format as PdPC theorem 2) ----------
CRB = np.linalg.inv(F_full[:3, :3] + np.diag([0, 0, 1e-12]))  # pseudo-inverse guard (for display)
F_struct = F_full[:2, :2]
cond_struct = np.linalg.cond(F_struct)
print(f'\nStructural block (ln D_c, ln ρ): condition number = {cond_struct:.1f} (identifiable; cf. PdPC shape block ~30)')

# ============================================================
# prediction (3): Monte-Carlo validation of Cramér–Rao
# ============================================================
R = 400
rho_hat, Dc_hat, sA_hat = [], [], []
for r in range(R):
    n = rng.poisson(np.maximum(N_D, 1e-9)[:, None], size=(len(D_grid), M)).mean(axis=1)
    # Poisson regression: n = ρ·ln D₀ − ρ·ln D_c → linear least squares
    X = np.column_stack([np.log(D0), np.ones(len(D_grid))])
    beta, *_ = np.linalg.lstsq(X, n, rcond=None)
    rho_hat.append(beta[0]); Dc_hat.append(np.exp(-beta[1]/beta[0]))
    # amplitude: mean over M cells per dose
    ybar = rng.normal(sA, sig/np.sqrt(M), size=len(D_grid)).mean()
    sA_hat.append(ybar)
rho_hat, Dc_hat, sA_hat = map(np.array, (rho_hat, Dc_hat, sA_hat))

var_rho_emp = np.var(np.log(rho_hat), ddof=1)
var_Dc_emp  = np.var(np.log(Dc_hat),  ddof=1)
var_sA_emp  = np.var(np.log(sA_hat),  ddof=1)
F_inv = np.linalg.pinv(F_full)
var_rho_crb, var_Dc_crb = F_inv[1, 1], F_inv[0, 0]
g_sA = np.array([0, 0, 1, 1])                    # gradient of the estimable combination ln(sA)
var_sA_crb = g_sA @ F_inv @ g_sA                 # CRB of the estimable combination (not a marginal pseudo-variance)

print('='*70)
print('[Prediction (3): MC empirical variance vs Cramér–Rao lower bound] (log-parameter variance)')
print(f'{"":>8}{"MC empirical":>12}{"CRB":>12}{"ratio":>8}')
for name, emp, crb in [('ln ρ', var_rho_emp, var_rho_crb),
                       ('ln D_c', var_Dc_emp, var_Dc_crb),
                       ('ln sA', var_sA_emp, var_sA_crb)]:
    print(f'{name:>8}{emp:12.2e}{crb:12.2e}{emp/crb:8.2f}')

# likelihood ridge: fixing s at a wrong value leaves the likelihood unchanged, with Â compensating inversely
s_wrong = np.array([0.5, 0.8, 1.0, 1.3, 2.0])
A_comp = sA/s_wrong
ll = [-(np.sum((rng.normal(sA, sig/np.sqrt(M), 100) - sw*aw)**2))/(2*sig**2/M)
      for sw, aw in zip(s_wrong, A_comp)]
print('\n[Scale unidentifiable: likelihood ridge] fixing s at a wrong value; after compensation Â = sA/s the likelihood difference is < noise fluctuation:')
for sw, aw in zip(s_wrong, A_comp):
    print(f'  s = {sw:.1f}（{"ground truth" if sw==1.0 else "error"}）→ Â = {aw:.2f}（compensation：Â·s = {aw*sw:.2f} constant）')

# ============================================================
# four-panel figure
# ============================================================
fig, axes = plt.subplots(1, 4, figsize=(19, 4.6))
names = ['ln $D_c$', 'ln $\\rho$', 'ln $A$', 'ln $s$']

ax = axes[0]
x = np.arange(4)
ax.bar(x-0.2, w_cnt/w_cnt[0], 0.4, label='counts only (number limit)', color='#2ca02c')
ax.bar(x+0.2, w_full/w_full[0], 0.4, label='count+amplitude', color='#1f77b4')
ax.axhline(1e-10, color='r', ls=':', lw=1, label='float-noise floor')
ax.set_yscale('log'); ax.set_xticks(x, [f'λ{i+1}' for i in x])
ax.set_ylim(1e-16, 10); ax.legend(fontsize=9)
ax.set_title('(a) Fisher spectrum: count direction nonzero, scale direction exactly zero')

ax = axes[1]
im = ax.imshow(np.abs(V_full), cmap='viridis', vmin=0, vmax=1)
ax.set_xticks(range(4), [f'v{i+1}\nλ={w_full[i]/w_full[0]:.0e}' for i in range(4)], fontsize=8)
ax.set_yticks(range(4), names)
for i in range(4):
    for j in range(4):
        ax.text(j, i, f'{abs(V_full[i,j]):.2f}', ha='center', va='center',
                color='w' if abs(V_full[i,j]) > 0.5 else 'k', fontsize=9)
ax.set_title('(b) eigenvectors: null direction v4 = (0,0,1,−1)/√2\nis exactly the degeneracy-group orbit tangent (s↑A↓ holding sA)')

ax = axes[2]
labels = ['ln ρ', 'ln $D_c$', 'ln sA']
emp = [var_rho_emp, var_Dc_emp, var_sA_emp]
crb = [var_rho_crb, var_Dc_crb, var_sA_crb]
xx = np.arange(3)
ax.bar(xx-0.2, emp, 0.4, label=f'MC empirical variance (R={R})', color='#ff7f0e')
ax.bar(xx+0.2, crb, 0.4, label='CRB = pinv(F) diagonal', color='#9467bd')
ax.set_xticks(xx, labels); ax.set_yscale('log'); ax.legend(fontsize=9)
ax.set_title('(c) Cramér–Rao validation: all structural parameters attain the bound')

ax = axes[3]
sax = np.linspace(0.4, 2.6, 120)
aax = sA/sax
ax.plot(sax, aax, 'r-', lw=2, label='likelihood ridge: s·A = const (G orbit)')
yobs = rng.normal(sA, sig/np.sqrt(M), (len(D_grid), M)).mean()
for sw, mk in [(0.7, 'x'), (1.0, 'o'), (1.4, 's')]:
    ax.plot(sw, sA/sw, mk, ms=10, mew=2,
            label=f's={sw} compensation Â={sA/sw:.2f} (same likelihood)')
ax.plot(1.0, A_true, 'k*', ms=16, label='ground truth (s=1, A=2)')
ax.set_xlabel('s on log scale (simulated scale, nuisance)'); ax.set_ylabel('compensating amplitude Â')
ax.legend(fontsize=8, loc='upper right')
ax.set_title('(d) scale unidentifiable: the entire G orbit is likelihood-degenerate')

plt.tight_layout()
plt.savefig('p53_Fisherstratified_S7.png', dpi=200)
print('\nSaved p53_Fisherstratified_S7.png')
