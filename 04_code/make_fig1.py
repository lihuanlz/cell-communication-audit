# -*- coding: utf-8 -*-
"""Fig. 1 | One audit, four systems, one ledger. (Nature-style, SVG + PNG twin)
All numbers transcribed from Nature_SI_v05_2026-09-23.md (S2, S7.0) and
Nature_main_v05_2026-09-23.md. ASCII-only labels, no em-dashes.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
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

# Okabe-Ito colorblind-safe palette
C_HIT  = '#009E73'   # bluish green
C_INT  = '#E69F00'   # orange/amber
C_FALS = '#D55E00'   # vermillion
C_LIM  = '#999999'   # grey
C_BLUE = '#0072B2'

OUTDIR = os.path.dirname(os.path.abspath(__file__))
MM = 1/25.4

def solve_u(xi, kappa):
    xi = np.atleast_1d(np.asarray(xi, float))
    f = lambda u, x: x*(1-u)*(kappa+u) - u*(kappa+1-u)
    out = np.empty_like(xi)
    for i, x in enumerate(xi.flat):
        out.flat[i] = brentq(f, 1e-16, 1-1e-16, args=(x,), xtol=1e-14, rtol=1e-14)
    return out

fig = plt.figure(figsize=(183*MM, 235*MM))

# ============================================================
# Panel a: five-step audit discipline + degeneracy-group inset
# ============================================================
axa = fig.add_axes([0.015, 0.745, 0.97, 0.235]); axa.axis('off')
axa.set_xlim(0, 1); axa.set_ylim(0, 1)

steps = [
    ("1. Derive the\ndegeneracy group", "orbit invariants =\nreportable content"),
    ("2. Identify invariant\nstatistics", "counts, ratios,\nordinals, event times"),
    ("3. Pre-register\ndecision lines", "frozen in writing\nbefore unblinding"),
    ("4. Blind scoring,\nfrozen seed", "SEED 20260814 (P2/P3),\n20260815 (P4-P7);\nbootstrap x2000"),
    ("5. Double-recorded\nverdict / erratum", "never-move rule;\n44b, 47b, 49b"),
]
n = len(steps)
bw, bh = 0.115, 0.34
ytop = 0.97
xs = np.linspace(0.005, 1-bw-0.005, n)
for i, ((title, sub), x) in enumerate(zip(steps, xs)):
    box = FancyBboxPatch((x, ytop-bh), bw, bh,
                         boxstyle="round,pad=0.008,rounding_size=0.012",
                         fc='#EAF1F8', ec=C_BLUE, lw=0.8)
    axa.add_patch(box)
    axa.text(x+bw/2, ytop-bh/2+0.065, title, ha='center', va='center',
             fontsize=6.3, fontweight='bold', color='#113a5d')
    axa.text(x+bw/2, ytop-bh/2-0.085, sub, ha='center', va='center',
             fontsize=5.0, color='#333333')
    if i < n-1:
        axa.add_patch(FancyArrowPatch((x+bw+0.004, ytop-bh/2),
                                      (xs[i+1]-0.004, ytop-bh/2),
                                      arrowstyle='-|>', mutation_scale=7,
                                      color=C_BLUE, lw=0.9))
axa.text(0.005, 0.52, "one discipline, applied identically to every system:",
         fontsize=6, style='italic', color='#444444', va='center')

# --- inset: degeneracy-group concept ---
iy0, ih = 0.03, 0.40
px0, pw = 0.015, 0.28
axa.add_patch(Rectangle((px0, iy0+0.03), pw, ih-0.10, fc='#FBFBFB', ec='0.55', lw=0.6))
axa.text(px0+pw/2, iy0+ih-0.035, 'parameter space (6 physical parameters)',
         ha='center', va='center', fontsize=5.6, color='#333333')
orbit_pts = [(px0+0.045, iy0+0.09), (px0+0.105, iy0+0.145),
             (px0+0.165, iy0+0.20), (px0+0.225, iy0+0.255)]
for j, (qx, qy) in enumerate(orbit_pts):
    axa.plot([qx], [qy], 'o', ms=4.5, mfc='white', mec=C_FALS, mew=1.1, zorder=5)
    axa.text(qx+0.008, qy+0.022, f'$\\theta_{j+1}$', fontsize=5.4, color=C_FALS)
axa.plot([orbit_pts[0][0]-0.012, orbit_pts[-1][0]+0.025],
         [orbit_pts[0][1]-0.010, orbit_pts[-1][1]+0.020],
         ls=':', lw=0.7, color='0.6', zorder=4)
axa.text(px0+pw/2, iy0+0.048, 'one group orbit: same $(\\xi,\\ \\kappa)$',
         ha='center', fontsize=5.2, color='0.35')

cx0, cw = 0.64, 0.34
axa.add_patch(Rectangle((cx0, iy0+0.03), cw, ih-0.10, fc='white', ec='0.55', lw=0.6))
axa.text(cx0+cw/2, iy0+ih-0.035, 'steady-state output curve',
         ha='center', va='center', fontsize=5.6, color='#333333')
xx = np.linspace(0, 1, 300)
uu = solve_u(np.exp(xx*6-3), 0.05)
axa.plot(cx0+0.045+xx*(cw-0.11), iy0+0.09+uu*(ih-0.26), color=C_BLUE, lw=1.2)
for (qx, qy) in orbit_pts:
    axa.add_patch(FancyArrowPatch((qx+0.025, qy), (cx0+0.015, iy0+0.16),
                                  arrowstyle='-|>', mutation_scale=5.5,
                                  color='0.45', lw=0.7))
axa.text((px0+pw+cx0)/2-0.01, iy0+0.26, 'group action', ha='center',
         fontsize=5.2, style='italic', color='0.4')
axa.text(cx0+cw/2, iy0+0.052, 'four parameter sets, one curve',
         ha='center', va='bottom', fontsize=5.4, color=C_FALS)
axa.text(cx0+cw/2, iy0-0.015, 'max $|\\Delta u|$ = 0 ($K_m,S_T$ x5); '
         '$1.7\\times10^{-18}$ ($k,E$ rescaling): machine precision',
         ha='center', va='top', fontsize=5.2, color=C_FALS)

axa.text(-0.012, 1.0, 'a', transform=axa.transAxes, fontsize=10,
         fontweight='bold', va='top')

# ============================================================
# Panel b: synthesis matrix (4 systems x 6 criteria)
# ============================================================
axb = fig.add_axes([0.015, 0.415, 0.97, 0.265]); axb.axis('off')
axb.set_xlim(0, 1); axb.set_ylim(0, 1)

rows = ['PdPC\nstatic sensor', 'p53\npulse encoder',
        'chemotaxis\nreceptor layer', 'GPCR\nfield meta-audit']
groups = [('Absolute scale', '#F4E3E3'), ('Scale-free', '#E3F0E8')]
cols = ['identifiable', 'internally\nconsistent', 'information-\nbearing',
        'identifiable', 'consistent', 'information-\nbearing']

NO, YES = 'NO', 'YES'
cells = [
    [(C_FALS, NO,  '4D degeneracy group:\n6 params -> 2 combos'),
     (C_FALS, NO,  'scale block: exact zero\neigenvalue ($\\chi\\sim10^{300}$)'),
     (C_FALS, NO,  'no information on\nabsolute scales'),
     (C_HIT,  YES, '$\\kappa$ exact from shape:\n$n_H = 1 + 1/(2\\kappa)$'),
     (C_HIT,  YES, 'closed form to machine\nprecision'),
     (C_HIT,  YES, 'interior optimum $\\kappa^*$:\n2.85-3.04 bits')],
    [(C_FALS, NO,  'amplitude: 5.5% of\nparameter information'),
     (C_FALS, NO,  'gain dispersion: CV(A)\n0.002 -> 0.639'),
     (C_FALS, NO,  'two-layer death of the\namplitude channel'),
     (C_HIT,  YES, 'timing channels carry\n94.2% of information'),
     (C_HIT,  YES, 'CV(T) pinned at 0.016\nunder gain dispersion'),
     (C_HIT,  YES, '$I(D;N)$ = 3.07 of 3.17\nbits (96.8%)')],
    [(C_FALS, NO,  'amplitude and midpoint\nfits irreconcilable'),
     (C_FALS, NO,  'Pareto rupture: no fit with\n$R^2$>0.9 and <0.1 dex'),
     (C_FALS, NO,  'blind amplitude predictions\nmiss (0.105, 0.126 vs 0.08)'),
     (C_HIT,  YES, 'per-cell midpoint\n$K_{1/2}(B)$'),
     (C_HIT,  YES, 'Weber line exact:\nmax dev. 0.10 dex'),
     (C_HIT,  YES, 'blind midpoint predictions\nhit (+0.014 to +0.038 dex)')],
    [(C_FALS, NO,  'functional vs binding:\nup to -3.5 dex (21-23$\\sigma$)'),
     (C_FALS, NO,  '21.5% of 17,987 pairs fail\nidentities >0.3 dex'),
     (C_FALS, NO,  'fractures at load-bearing\npositions'),
     (C_HIT,  YES, 'identity-based checks\n(scale-free)'),
     (C_HIT,  YES, 'pipeline exact:\n3,931/3,931 closures'),
     (C_HIT,  YES, 'repair named by model:\n$\\Delta$logML = +128.2')],
]

tbl_x0, tbl_x1 = 0.115, 0.995
tbl_y0, tbl_y1 = 0.06, 0.76
ncol, nrow = 6, 4
cw_ = (tbl_x1-tbl_x0)/ncol
rh = (tbl_y1-tbl_y0)/nrow
for gi, (gname, gtint) in enumerate(groups):
    gx0 = tbl_x0 + gi*3*cw_
    axb.add_patch(Rectangle((gx0+0.004, tbl_y1+0.095), 3*cw_-0.008, 0.085,
                            fc=gtint, ec='0.4', lw=0.6))
    axb.text(gx0+1.5*cw_, tbl_y1+0.137, gname, ha='center', va='center',
             fontsize=6.5, fontweight='bold')
    for j in range(3):
        axb.text(gx0+(j+0.5)*cw_, tbl_y1+0.048, cols[gi*3+j], ha='center',
                 va='center', fontsize=5.2)
for i, rname in enumerate(rows):
    y1 = tbl_y1 - i*rh
    axb.text(tbl_x0-0.012, y1-rh/2, rname, ha='right', va='center',
             fontsize=6.0, fontweight='bold')
    for j in range(ncol):
        col, word, txt = cells[i][j]
        x0 = tbl_x0 + j*cw_
        axb.add_patch(Rectangle((x0+0.005, y1-rh+0.008), cw_-0.010, rh-0.016,
                                fc=col, ec='none', alpha=0.16))
        axb.add_patch(Rectangle((x0+0.005, y1-rh+0.008), cw_-0.010, rh-0.016,
                                fc='none', ec=col, lw=0.8))
        axb.text(x0+cw_/2, y1-rh*0.24, word, ha='center', va='center',
                 fontsize=6.2, fontweight='bold', color=col)
        axb.text(x0+cw_/2, y1-rh*0.62, txt, ha='center', va='center',
                 fontsize=4.9, color='#222222')
axb.text(tbl_x0+1.5*cw_, tbl_y0-0.030,
         'unidentifiable / inconsistent / unused', ha='center', va='top',
         fontsize=5.2, style='italic', color=C_FALS)
axb.text(tbl_x0+4.5*cw_, tbl_y0-0.030,
         'identifiable / exact / information-bearing', ha='center', va='top',
         fontsize=5.2, style='italic', color=C_HIT)
axb.text(-0.012, 1.06, 'b', transform=axb.transAxes, fontsize=10,
         fontweight='bold', va='top')

# ============================================================
# Panel c: blind-adjudication ledger (14 clauses)
# ============================================================
axc = fig.add_axes([0.175, 0.075, 0.475, 0.285])

M = lambda v, line: v/line - 1.0
rows_c = [
 ('P2-1 Msn2: size-tercile penalty (mol)', C_FALS,
  [(M(0.0294, 0.03), 'o')], None,
  '$\\Delta_{mol}$ max +0.0294 vs line +0.03;\n0/7 promoters pass'),
 ('P2-2 Msn2: AU-channel penalty bounded', C_FALS,
  [(-1.05, 'D')], None,
  '0/7 promoters pass (line: 6/7);\nposition conventional'),
 ('P2-3 Msn2: event-time stratification-free', C_HIT,
  [(-M(0.0102, 0.01), 'o')], None,
  '6/7 pass; worst $|\\Delta|$ = 0.0102 vs line\n0.01 (7/7 at k=2, 6/7 at k=4)'),
 ('P3-1 NF-$\\kappa$B grad: dose ordinality', C_INT,
  [(M(0.6548, 0.60), 'o'), (M(0.5184, 0.60), 'o')], None,
  'AUC 0.655 (pass) / 0.518 (fail)\nvs line 0.60'),
 ('P3-2 NF-$\\kappa$B grad: spatial ordinality', C_HIT,
  [(M(0.9070, 0.55), 'o'), (M(0.8193, 0.55), 'o')], None,
  'AUC 0.907 / 0.819 vs line 0.55,\nboth pass'),
 ('P3-3 NF-$\\kappa$B grad: duration axis', C_FALS,
  [(M(0.5911, 0.60), 'o'), (M(0.5214, 0.60), 'o')], None,
  'AUC 0.591 / 0.521, both below\nline 0.60'),
 ('P4-1 ERK-KTR: event-time decoding', C_FALS,
  [(M(v, 0.60), 'o') for v in (0.5819, 0.5633, 0.5636, 0.4695,
                                0.5561, 0.5663, 0.5377, 0.5096)], None,
  'all 8 AUCs below 0.60\n(max 0.582; one 0.470 < 0.5)'),
 ('P4-2 ERK-KTR: saturation ordering', C_INT,
  [(M(0.0468, 0.05), 'o')], (M(0.032, 0.05), M(0.062, 0.05)),
  'D = +0.0468 [+0.032, +0.062]\nvs hit line 0.05'),
 ('P4-3 ERK-KTR: size-gain invariance', C_FALS,
  [((4-8)/8, 'o')], None,
  '5/10 units $|\\Delta|\\geq$0.03 (falsif. line\n3/10); hit line 8/10, reached 4/10'),
 ('P5-1 GPCR-Ca: dead-zone clause', C_FALS,
  [(-M(0.629, 0.60), 'o'), (-M(0.659, 0.60), 'o'),
   ((0.60-0.568)/0.60, 'o'), ((0.60-0.520)/0.60, 'o')], None,
  '2/4 units out of band (0.629, 0.659);\nband [0.40, 0.60)'),
 ('P5-2 GPCR-Ca: live-zone counting', C_INT,
  [(M(v, 0.60), 'o') for v in (0.629, 0.659, 0.570, 0.522)], None,
  '2/4 units at line (0.629, 0.659;\nbelow: 0.570, 0.522)'),
 ('P5-3 GPCR-Ca: paid-peak control', C_LIM,
  [(M(v, 0.60), 'o') for v in (0.731, 0.744, 0.551, 0.523)], None,
  '2/4 units >= 0.60 (0.731, 0.744);\nneeded 3/4'),
 ('P6 NF-$\\kappa$B seq: fold vs absolute', C_INT,
  [(M(-0.0014, 0.15), 'o')], (M(-0.2518, 0.15), M(0.1533, 0.15)),
  '$\\Delta R^2$ = -0.0014 [-0.252, +0.153],\n18 units; CI crosses zero'),
 ('P7 chemo ladder: fold/Weber vs total', C_INT,
  [(M(-0.0909, 0.15), 'o'), (M(0.1181, 0.15), 's')],
  (M(-0.4446, 0.15), M(0.2766, 0.15)),
  'erratum $\\Delta R^2$ = -0.091 [-0.445, +0.277],\n30 units; first pass +0.118 (square)'),
]

XL = 1.30
nrow = len(rows_c)
axc.set_ylim(-0.8, nrow-0.4)
axc.set_xlim(-XL, XL)
axc.axvspan(-XL, 0, color=C_FALS, alpha=0.05)
axc.axvspan(0, XL, color=C_HIT, alpha=0.05)
axc.axvline(0, color='k', lw=0.9)

tickpos, ticklab = [], []
gaps = {3: 'P2', 6: 'P3', 9: 'P4', 12: 'P5', 13: 'P6'}  # row index -> end of group
for i, (label, col, pts, ci, note) in enumerate(rows_c):
    y = nrow-1-i
    tickpos.append(y); ticklab.append(label)
    if ci is not None:
        lo, hi = ci
        xa, xb = max(lo, -XL), min(hi, XL)
        axc.plot([xa, xb], [y, y], color=col, lw=0.9, zorder=3)
        if lo < -XL:
            axc.annotate('', xy=(-XL+0.02, y), xytext=(-XL+0.16, y),
                         arrowprops=dict(arrowstyle='-|>', color=col, lw=0.9))
        if hi > XL:
            axc.annotate('', xy=(XL-0.02, y), xytext=(XL-0.16, y),
                         arrowprops=dict(arrowstyle='-|>', color=col, lw=0.9))
    for x, mk in pts:
        ms = 4.2 if mk == 'o' else 4.6
        axc.plot([x], [y], marker=mk, ms=ms, mfc=col, mec='k', mew=0.4,
                 zorder=5, linestyle='none')
    axc.text(XL+0.07, y, note, fontsize=4.9, va='center', color='#222222',
             clip_on=False)
    if i in (2, 5, 8, 11, 12):
        axc.axhline(y-0.5, color='0.78', lw=0.5)

axc.set_yticks(tickpos); axc.set_yticklabels(ticklab, fontsize=5.4)
axc.set_xticks([-1, -0.5, 0, 0.5, 1])
axc.set_xticklabels(['-1', '-0.5', 'frozen line', '+0.5', '+1'], fontsize=5.6)
axc.set_xlabel('signed distance to frozen decision line, in units of the line (+ = hit side)',
               fontsize=6.3, labelpad=2)
axc.tick_params(length=2)
axc.xaxis.set_label_coords(0.5, -0.13)
for s in ('top', 'right'):
    axc.spines[s].set_visible(False)

handles = [plt.Line2D([], [], marker='o', ls='none', mfc=c, mec='k', mew=0.4, ms=4.5)
           for c in (C_HIT, C_INT, C_FALS, C_LIM)]
fig.legend(handles, ['hit (2)', 'intermediate (5)', 'falsified (6)',
                     'dataset-level limitation (1)'],
           loc='lower left', bbox_to_anchor=(0.655, 0.362), ncol=2, fontsize=5.4,
           frameon=False, handletextpad=0.15, columnspacing=0.8)
axc.set_title('Totals: 2 hits / 5 intermediates / 6 falsifications / 1 dataset-level limitation\n'
              '(14 pre-registered clauses, 7 adjudications P2-P7, 6 datasets)',
              fontsize=6.0, fontweight='bold', loc='left', pad=14)
axc.text(-0.375, 1.085, 'c', transform=axc.transAxes, fontsize=10,
         fontweight='bold', va='top')

fig.savefig(os.path.join(OUTDIR, 'Fig1_ledger.svg'))
fig.savefig(os.path.join(OUTDIR, 'Fig1_ledger.png'), dpi=200)
print('Fig1 written.')
