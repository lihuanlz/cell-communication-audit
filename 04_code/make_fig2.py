# -*- coding: utf-8 -*-
"""Fig. 2 | The static sensor (PdPC). (Nature-style, SVG + PNG twin)
Numbers transcribed from _archive_旧版本/Nature_SI_v04_2026-09-23.md S2 (tables in S2.3, S2.6)
and verified against code 4 functional forms. ASCII-only labels, no em-dashes.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7,
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.labelsize': 6, 'ytick.labelsize': 6,
    'svg.fonttype': 'none',
    'mathtext.fontset': 'dejavusans',
})

C_BLUE, C_VERM, C_GRN, C_ORG, C_SKY, C_GREY = ('#0072B2', '#D55E00', '#009E73',
                                               '#E69F00', '#56B4E9', '#999999')
OUTDIR = os.path.dirname(os.path.abspath(__file__))
MM = 1/25.4
SIG = 0.03

def solve_u(xi, kappa):
    xi = np.atleast_1d(np.asarray(xi, float))
    if kappa < 1e-12:
        return np.where(xi < 1, 0.0, np.where(xi > 1, 1.0, 0.5))
    f = lambda u, x: x*(1-u)*(kappa+u) - u*(kappa+1-u)
    out = np.empty_like(xi)
    for i, x in enumerate(xi.flat):
        out.flat[i] = brentq(f, 1e-16, 1-1e-16, args=(x,), xtol=1e-14, rtol=1e-14)
    return out

def dxi_du(u, kappa):
    xi = u/(1-u)*(kappa+1-u)/(kappa+u)
    return xi*(1/u - 1/(kappa+1-u) + 1/(1-u) - 1/(kappa+u))

def du_dk(u, xi, kappa):
    Gk = xi*(1-u) - u
    Gu = xi*(1-kappa-2*u) - (kappa+1-2*u)
    return -Gk/Gu

fig = plt.figure(figsize=(183*MM, 118*MM))
gs = fig.add_gridspec(2, 6, top=0.86, bottom=0.135, left=0.065, right=0.985,
                      hspace=0.85, wspace=0.9)

def panel_letter(ax, s, dx=-0.16):
    ax.text(dx, 1.04, s, transform=ax.transAxes, fontsize=10,
            fontweight='bold', va='top')

# ---------------- Panel a: steady-state curves across kappa ----------------
axa = fig.add_subplot(gs[0, 0:2])
xi_f = np.logspace(-3, 3, 600)
kaps_a = [0.001, 0.01, 0.1, 1.0]
cols_a = [C_VERM, C_ORG, C_GRN, C_BLUE]
for k, c in zip(kaps_a, cols_a):
    nH = 1 + 1/(2*k)
    nH_s = f'{nH:.4g}'
    axa.semilogx(xi_f, solve_u(xi_f, k), lw=1.4, color=c,
                 label=f'$\\kappa$ = {k:g}  ($n_H$ = {nH_s})')
axa.set_xlabel('input activity ratio $\\xi$')
axa.set_ylabel('modified fraction $u^*$')
axa.set_ylim(0, 1.02)
axa.legend(fontsize=5.4, loc='lower right', frameon=False, handlelength=1.4)
axa.set_title('symmetric Michaelis steady state', fontsize=6.3)
panel_letter(axa, 'a')

# ---------------- Panel b: degeneracy-group invariance ----------------
axb = fig.add_subplot(gs[0, 2:4])
base = dict(k1=1.0, E1T=1.0, k2=1.0, E2T=1.0, Km=0.1, ST=1.0)
x_grid = np.logspace(-2, 2, 200)
def tcurve(p):
    xi = p['k1']*p['E1T']*x_grid/(p['k2']*p['E2T'])
    return solve_u(xi, p['Km']/p['ST'])
sets = [
    ('$\\theta_1$ (1,1,1,1,0.1,1)', base, dict(color='k', lw=2.6, ls='-', marker=None)),
    ('$\\theta_2$ $k_1,k_2$ x2', {**base, 'k1':2.0, 'k2':2.0},
     dict(color=C_VERM, lw=1.1, ls='--', marker='o')),
    ('$\\theta_3$ $K_m,S_T$ x5', {**base, 'Km':0.5, 'ST':5.0},
     dict(color=C_BLUE, lw=1.1, ls='-.', marker='s')),
    ('$\\theta_4$ $k_1,E_{2T}$ x3; $K_m,S_T$ x2', {**base, 'k1':3.0, 'E2T':3.0, 'Km':0.2, 'ST':2.0},
     dict(color=C_GRN, lw=1.1, ls=':', marker='^')),
]
for name, p, st in sets:
    u = tcurve(p)
    if st['marker'] is None:
        axb.semilogx(x_grid, u, color=st['color'], lw=st['lw'], ls=st['ls'], label=name)
    else:
        axb.semilogx(x_grid[::14], u[::14], color=st['color'], lw=st['lw'],
                     ls=st['ls'], marker=st['marker'], ms=4, mfc='none',
                     label=name)
axb.set_xlabel('titration variable $x$')
axb.set_ylabel('modified fraction $u^*$')
axb.set_ylim(0, 1.02)
axb.legend(fontsize=5.0, loc='upper left', frameon=False, handlelength=1.6)
axb.annotate('four parameter sets, one orbit, one curve\n'
             'max $|\\Delta u|$ = 0 ($K_m,S_T$ x5)\n'
             'max $|\\Delta u|$ = $1.7\\times10^{-18}$ ($k,E$)',
             xy=(0.98, 0.03), xycoords='axes fraction', ha='right', va='bottom',
             fontsize=5.2, color=C_VERM)
axb.set_title('degeneracy-group invariance (200-pt grid)', fontsize=6.3)
panel_letter(axb, 'b')

# ---------------- Panel c: Fisher eigenspectrum ----------------
axc = fig.add_subplot(gs[0, 4:6])
LNXI_D = np.linspace(-1.5*np.log(10), 1.5*np.log(10), 20)
XI_D = np.exp(LNXI_D)
kap_c = 0.05
u_d = solve_u(XI_D, kap_c)
a_ = XI_D/dxi_du(u_d, kap_c)
b_ = kap_c*du_dk(u_d, XI_D, kap_c)
Fs = np.column_stack([a_, b_]); Fs = Fs.T @ Fs / SIG**2
Fa = np.column_stack([a_, b_, -b_]); Fa = Fa.T @ Fa / SIG**2
ev_s = np.linalg.eigvalsh(Fs)[::-1]
ev_a = np.linalg.eigvalsh(Fa)[::-1]
print('shape eigenvalues:', ev_s, ' chi_shape =', ev_s[0]/ev_s[1])
print('abs eigenvalues:', ev_a, ' chi_abs(num) =', ev_a[0]/abs(ev_a[2]))
labels = ['$\\lambda_1$', '$\\lambda_2$', '$\\lambda_1$', '$\\lambda_2$',
          '$\\lambda_3$']
vals = [ev_s[0], ev_s[1], ev_a[0], ev_a[1], max(abs(ev_a[2]), 1e-14)]
colors = [C_BLUE, C_BLUE, C_VERM, C_VERM, C_VERM]
xp = np.arange(5)
axc.bar(xp, np.log10(vals), width=0.62, color=colors, alpha=0.85)
axc.set_ylabel('$\\log_{10}$ Fisher eigenvalue')
axc.set_xticks(xp); axc.set_xticklabels(labels, fontsize=6)
axc.set_ylim(-23, 7.4)
axc.axvspan(-0.5, 1.5, color=C_BLUE, alpha=0.06)
axc.axvspan(1.5, 4.5, color=C_VERM, alpha=0.06)
axc.text(0.5, 7.2, 'shape tier\n($\\ln c, \\ln\\kappa$)',
         ha='center', va='top', fontsize=5.4, color=C_BLUE, fontweight='bold')
axc.text(3, 7.2, 'absolute tier\n($\\ln c, \\ln K_m, \\ln S_T$)',
         ha='center', va='top', fontsize=5.4, color=C_VERM, fontweight='bold')
axc.text(0.5, -2.2, '$\\chi \\approx$ 19-34', ha='center', va='top',
         fontsize=5.6, color=C_BLUE, fontweight='bold')
axc.text(2.72, -2.2, 'exact zero\neigenvalue along orbit', ha='center', va='top',
         fontsize=5.2, color=C_VERM, fontweight='bold')
axc.annotate('', xy=(4, -21.5), xytext=(4, -13.2),
             arrowprops=dict(arrowstyle='-|>', color=C_VERM, lw=1.0))
axc.text(3.55, -16.0,
         '$\\lambda_{min}$ at float floor ~$10^{-14}$:\n'
         'effective $\\chi$ saturates at $10^{17}$-$10^{19}$;\n'
         'extrapolated $\\chi\\sim10^{300}$ (deep saturation)',
         fontsize=5.0, va='center', ha='right', color=C_VERM)
axc.set_title('Fisher eigenspectrum at $\\kappa$ = 0.05 (schematic tiers)', fontsize=6.3)
panel_letter(axc, 'c')

# ---------------- Panel d: Hill-slope closed form ----------------
axd = fig.add_subplot(gs[1, 0:3])
kk = np.logspace(-3.2, 0.2, 200)
axd.loglog(kk, 1+1/(2*kk), 'k-', lw=1.4, label='closed form  $n_H = 1 + 1/(2\\kappa)$')
kap_tab = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0]
nH_tab  = [501, 101, 51, 11, 6, 2, 1.5]   # SI S2.3 table (= closed form)
axd.loglog(kap_tab, nH_tab, 'o', ms=5, mfc='none', mec=C_VERM, mew=1.2,
           ls='none', label='numerical logit slope (archived)')
axd.axvspan(3e-4, 0.01, color=C_ORG, alpha=0.12)
axd.text(3.5e-3, 18, 'one-sided-bound regime:\nreport $\\kappa \\leq \\kappa_{max}$,\nnever a point estimate',
         fontsize=5.4, ha='center', color='#8a5a00')
axd.annotate('$\\kappa = 10^{-3}$: $n_H$ = 501,\n95% CI inflation x2920',
             xy=(0.0015, 380), xytext=(0.02, 120),
             fontsize=5.6, color=C_VERM,
             arrowprops=dict(arrowstyle='-|>', color=C_VERM, lw=0.8))
axd.set_xlabel('$\\kappa = K_m/S_T$')
axd.set_ylabel('Hill slope $n_H$')
axd.legend(fontsize=5.6, loc='lower left', frameon=False)
axd.set_title('Hill-slope closed form vs numerics', fontsize=6.3)
panel_letter(axd, 'd', dx=-0.11)

# ---------------- Panel e: mutual information vs kappa ----------------
axe = fig.add_subplot(gs[1, 3:6])
kap_mi = [0.002, 0.01, 0.05, 0.1, 0.3, 1, 3, 10]
I_narrow = [1.42, 2.12, 2.84, 2.85, 2.44, 1.87, 1.58, 1.46]  # SI S2.6, sigma_xi=0.3
I_wide   = [1.11, 1.41, 2.07, 2.42, 2.85, 3.03, 3.04, 3.02]  # SI S2.6, sigma_xi=1.5
axe.axhspan(1.1, 1.4, color=C_GREY, alpha=0.18)
axe.text(1.55e-3, 0.78, 'zero-order limit ($\\kappa\\to0$): 1.1-1.4 bits;\n'
         'a decision element, not a measurement element',
         fontsize=5.2, va='center', color='0.30')
axe.semilogx(kap_mi, I_narrow, '^-', color=C_BLUE, ms=5, lw=1.2,
             label='narrow input, $\\sigma_\\xi$ = 0.3')
axe.semilogx(kap_mi, I_wide, 's-', color=C_VERM, ms=4.6, lw=1.2,
             label='wide input, $\\sigma_\\xi$ = 1.5')
axe.axvline(0.1, color=C_BLUE, ls=':', lw=0.9)
axe.axvline(3, color=C_VERM, ls=':', lw=0.9)
axe.annotate('$\\kappa^* \\approx$ 0.05-0.1\n2.85 bits', xy=(0.1, 2.85),
             xytext=(0.0035, 3.25), fontsize=5.4, color=C_BLUE,
             arrowprops=dict(arrowstyle='-|>', color=C_BLUE, lw=0.8))
axe.annotate('$\\kappa^* \\approx$ 1-3\n3.04 bits', xy=(3, 3.04),
             xytext=(0.45, 3.32), fontsize=5.4, color=C_VERM,
             arrowprops=dict(arrowstyle='-|>', color=C_VERM, lw=0.8))
axe.set_xlabel('$\\kappa = K_m/S_T$')
axe.set_ylabel('$I(\\xi; u)$  (bits)')
axe.set_ylim(0.55, 3.6)
axe.set_xlim(1.4e-3, 14)
axe.legend(fontsize=5.6, loc='lower right', frameon=False)
axe.set_title('interior information optimum; zero-order limit is poor', fontsize=6.3)
panel_letter(axe, 'e', dx=-0.11)

fig.savefig(os.path.join(OUTDIR, 'Fig2_static_sensor.svg'))
fig.savefig(os.path.join(OUTDIR, 'Fig2_static_sensor.png'), dpi=200)
print('Fig2 written.')
