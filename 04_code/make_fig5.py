# -*- coding: utf-8 -*-
"""
Fig. 5 | The receptor-layer rupture and the prospective test.
All numbers transcribed from _archive_旧版本/Nature_SI_v04_2026-09-23.md S8, code 54/62 docstrings,
and code70 run output. Nature-style: Arial/Helvetica fallback DejaVu Sans, 7pt,
svg.fonttype='none', axes linewidth 0.6, Okabe-Ito palette.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
import matplotlib.gridspec as gridspec

plt.rcParams.update({
    "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "svg.fonttype": "none",
    "axes.unicode_minus": False,
})

# Okabe-Ito
BLUE, VERM, GREEN, ORANGE, GREY, SKY = ("#0072B2", "#D55E00", "#009E73",
                                        "#E69F00", "#7F7F7F", "#56B4E9")

OUT = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911\05_主线纲领与设计\论文_细胞通讯审计_2026-09-23\figures_svg"

fig = plt.figure(figsize=(7.2, 7.4))
gs = gridspec.GridSpec(3, 2, height_ratios=[1.05, 1.25, 1.0],
                       hspace=0.42, wspace=0.30,
                       left=0.075, right=0.975, top=0.965, bottom=0.075)

def panel_label(ax, s, dx=-0.16, dy=1.04):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="top", ha="left")

# ---------------------------------------------------------------- panel a
# Pareto front of joint fits (code 54 / SI S8.2)
axa = fig.add_subplot(gs[0, 0])
# front points: (amplitude R2, K1/2 error dex, label)
front = [(0.980, 0.622, "N = 39"), (0.831, 0.22, ""), (0.511, 0.074, "N = 2.0")]
xs = [p[0] for p in front]; ys = [p[1] for p in front]
# smooth interpolating curve through the three documented points
tt = np.linspace(0, 1, 200)
# monotone-ish curve in log-lambda space: use quadratic in t through 3 pts
from numpy.polynomial import polynomial as P
cx = np.polyfit([0, 0.5, 1], xs, 2); cy = np.polyfit([0, 0.5, 1], ys, 2)
axa.plot(np.polyval(cx, tt), np.polyval(cy, tt), color=BLUE, lw=1.2, zorder=2)
# acceptance region: R2 > 0.9 and K1/2 error < 0.1
axa.add_patch(Rectangle((0.9, 0.0), 0.10, 0.10, facecolor=GREEN, alpha=0.18,
                        edgecolor=GREEN, lw=0.8, zorder=1))
axa.text(0.885, 0.05, "acceptance\nregion (empty)", fontsize=6, color=GREEN,
         ha="right", va="center")
axa.plot(xs, ys, "o", color=BLUE, ms=5, mfc="white", mew=1.2, zorder=3)
axa.annotate("N = 39\n(0.980, 0.622)", (0.980, 0.622), xytext=(0.80, 0.55),
             fontsize=6, color=BLUE,
             arrowprops=dict(arrowstyle="-", lw=0.5, color=BLUE))
axa.annotate("N = 2.0\n(0.511, 0.074)", (0.511, 0.074), xytext=(0.44, 0.24),
             fontsize=6, color=BLUE,
             arrowprops=dict(arrowstyle="-", lw=0.5, color=BLUE))
# best repair compromises (SI S8.3), all excluded
rep = [(0.953, 0.182, "M2c"), (0.950, 0.257, "V3"), (0.946, 0.179, "het.")]
for x, y, s in rep:
    axa.plot(x, y, "x", color=GREY, ms=4.5, mew=1.1, zorder=3)
axa.annotate("best repair\ncompromises\n(excluded)", (0.950, 0.257),
             xytext=(0.62, 0.40), fontsize=6, color=GREY, ha="left",
             arrowprops=dict(arrowstyle="-", lw=0.5, color=GREY))
axa.text(0.415, 0.66, "midpoint statistic\nrequires N ~ 2-4", fontsize=6,
         color="black", ha="left")
axa.text(0.80, 0.685, "amplitude statistic\nrequires N ~ 12-39", fontsize=6,
         color="black", ha="left")
axa.set_xlim(0.40, 1.03); axa.set_ylim(0.0, 0.75)
axa.set_xlabel("amplitude $R^2$")
axa.set_ylabel("$K_{1/2}$ error (dex)")
axa.set_title("Pareto front of joint fits", fontsize=7)
panel_label(axa, "a")

# ---------------------------------------------------------------- panel b
# Weber line (SI S8.4; code 54 K1/2 medians)
axb = fig.add_subplot(gs[0, 1])
B = np.array([0, 0.01, 0.1, 0.3, 1, 10, 100.0])
K = np.array([2.03, 2.90, 2.17, 2.58, 3.49, 13.92, 119.56])
Bx = np.where(B == 0, 0.006, B)          # plot B=0 at left edge
line = lambda b: 1.17 * (1.95 + b)
bx = np.logspace(np.log10(0.006), 2, 100)
axb.plot(bx, line(bx), color=VERM, lw=1.2, zorder=2,
         label="Weber line $K_{1/2} = 1.17\\,(1.95 + B)$ µM")
axb.plot(Bx, K, "o", color=BLUE, ms=4.5, mfc="white", mew=1.2, zorder=3,
         label="measured $K_{1/2}$ medians\n(57-268 cells per background)")
axb.set_xscale("log"); axb.set_yscale("log")
axb.set_xticks([0.006, 0.01, 0.1, 1, 10, 100])
axb.set_xticklabels(["0", "0.01", "0.1", "1", "10", "100"])
axb.set_yticks([1, 3, 10, 30, 100])
axb.set_yticklabels(["1", "3", "10", "30", "100"])
axb.set_xlim(0.005, 180); axb.set_ylim(1.2, 220)
axb.set_xlabel("background $B$ (µM MeAsp)")
axb.set_ylabel("$K_{1/2}$ median (µM)")
axb.legend(fontsize=5.6, loc="upper left", frameon=False, handlelength=1.4)
axb.set_title("Weber line across seven backgrounds", fontsize=7)
panel_label(axb, "b")
# inset residuals
axin = axb.inset_axes([0.40, 0.06, 0.57, 0.27])
res = np.log10(K / line(B))
axin.axhline(0, color="black", lw=0.5)
axin.plot(Bx, res, "o", color=BLUE, ms=3, mfc="white", mew=0.9)
axin.set_xscale("log")
axin.set_xlim(0.005, 180); axin.set_ylim(-0.14, 0.14)
axin.set_xticks([0.006, 0.1, 10]); axin.set_xticklabels(["0", "0.1", "10"], fontsize=5)
axin.set_yticks([-0.1, 0, 0.1]); axin.tick_params(labelsize=5, width=0.5, size=1.8)
axin.set_ylabel("residual (dex)", fontsize=5.5, labelpad=1)
for s in axin.spines.values():
    s.set_linewidth(0.5)
axin.text(0.03, 0.90, "max |resid| = 0.10 dex", fontsize=5.5,
          transform=axin.transAxes, ha="left", va="top")

# ---------------------------------------------------------------- panel c
# Eight repair paths, all excluded (SI S8.3)
axc = fig.add_subplot(gs[1, :])
axc.axis("off")
panel_label(axc, "c", dx=-0.02, dy=1.06)
repairs = [
    ("1. Two subpopulations\n(M2a / M2b)",
     "directional failure: alpha driven to 0.05-0.08,\nworse than single-species fit"),
    ("2. Imperfect adaptation\n(M2c, beta < 1)",
     "best 0.953 / 0.182 dex; gain from Ka to inf,\nnot beta; outside acceptance"),
    ("3. Protocol mismatch\n(step vs waveform)",
     "waveform FCD plateau only at B = 10-100 µM;\nexcluded as primary cause"),
    ("4. FRET readout\nnonlinearity (V1-V3)",
     "V3 picks q ~ 2.0 expansive; best 0.950 / 0.257;\nexcluded as sole cause"),
    ("5. TCS ligand depletion\n(kappa ~ 48)",
     "free-ligand solution leaves the front\ndigit-for-digit unchanged"),
    ("6. Background-dependent\ngain retuning",
     "K1/2 side 0.052 dex but R2 collapses to 0.68;\nnon-physical boundaries hit"),
    ("7. Heterogeneous population\n+ censoring (code 60)",
     "best 0.946 / 0.179; optimizer shrinks sigma_a\nto 0.11, opposite direction"),
    ("8. Asymmetric cooperativity\n(N ~ 2-4 vs 12-39)",
     "measured n(B) within single-N MWC envelope;\ntwo N values not required"),
]
for i, (title, reason) in enumerate(repairs):
    r, c = divmod(i, 4)
    x0 = 0.015 + c * 0.247; y0 = 0.54 - r * 0.52
    axc.add_patch(FancyBboxPatch((x0, y0), 0.225, 0.42,
                                 boxstyle="round,pad=0.008", fc="#F5F5F5",
                                 ec="0.4", lw=0.6, transform=axc.transAxes,
                                 clip_on=False))
    axc.text(x0 + 0.112, y0 + 0.345, title, transform=axc.transAxes,
             fontsize=6.0, fontweight="bold", ha="center", va="top")
    axc.text(x0 + 0.112, y0 + 0.175, reason, transform=axc.transAxes,
             fontsize=5.4, ha="center", va="top", color="0.25")
    axc.text(x0 + 0.112, y0 + 0.055, "EXCLUDED", transform=axc.transAxes,
             fontsize=7, fontweight="bold", color=VERM, ha="center", va="center",
             bbox=dict(boxstyle="round,pad=0.25", fc="none", ec=VERM, lw=0.9))
axc.set_title("Eight repair paths, all excluded (uniform criterion: $R^2$ > 0.9, "
              "$K_{1/2}$ error < 0.1 dex, $a_{max}$ ~ 1)", fontsize=7, pad=10)

# ---------------------------------------------------------------- panel d
# Same-source verdict (code 62; SI S8.5)
axd = fig.add_subplot(gs[2, 0])
ss = [(0.971, 0.573), (0.868, 0.205), (0.599, 0.071)]   # same-source
wp = [(0.976, 0.624), (0.848, 0.261), (0.525, 0.074)]   # whole-population control
axd.add_patch(Rectangle((0.9, 0.0), 0.10, 0.10, facecolor=GREEN, alpha=0.18,
                        edgecolor=GREEN, lw=0.8, zorder=1))
for pts, col, lab in [(wp, GREY, "whole-population control"),
                      (ss, BLUE, "same-source front (code 62)")]:
    x = [p[0] for p in pts]; y = [p[1] for p in pts]
    ccx = np.polyfit([0, 0.5, 1], x, 2); ccy = np.polyfit([0, 0.5, 1], y, 2)
    axd.plot(np.polyval(ccx, tt), np.polyval(ccy, tt), color=col, lw=1.1)
    axd.plot(x, y, "o", color=col, ms=4.5, mfc="white", mew=1.1, label=lab)
axd.text(0.948, 0.13, "acceptance\nregion", fontsize=6, color=GREEN,
         ha="center", va="bottom")
axd.text(0.44, 0.38, "front unmoved:\nrupture persists inside a\nstrictly single population",
         fontsize=6, color="black")
axd.annotate("lambda = 0 / 1 / 30", (0.868, 0.205), xytext=(0.60, 0.30),
             fontsize=6, color="0.25",
             arrowprops=dict(arrowstyle="-", lw=0.5, color="0.25"))
axd.set_xlim(0.40, 1.03); axd.set_ylim(0.0, 0.75)
axd.set_xlabel("amplitude $R^2$"); axd.set_ylabel("$K_{1/2}$ error (dex)")
axd.legend(fontsize=5.6, loc="upper left", frameon=False)
axd.set_title("Same-source verdict", fontsize=7)
panel_label(axd, "d")

# ---------------------------------------------------------------- panel e
# Prospective mixed-ligand arm (code 70; SI S8.6)
axe = fig.add_subplot(gs[2, 1])
# midpoint hits (frozen band +/-0.15 dex)
mid = [("V1", 0.014), ("V2", 0.038), ("V2_alt", 0.031)]
amp = [("V1", 0.105), ("V2", 0.126), ("V2_alt*", 0.069)]
axe.add_patch(Rectangle((0.55, 0), 2.9, 0.15, facecolor=GREEN, alpha=0.15,
                        edgecolor="none", zorder=0))
axe.add_patch(Rectangle((3.65, 0), 2.9, 0.08, facecolor=GREEN, alpha=0.15,
                        edgecolor="none", zorder=0))
axe.axhline(0.15, xmin=0.02, xmax=0.46, color=GREEN, lw=0.8, ls="--")
axe.axhline(0.08, xmin=0.54, xmax=0.98, color=VERM, lw=0.8, ls="--")
for i, (s, v) in enumerate(mid):
    axe.plot(1 + i, v, "o", color=GREEN, ms=6, mfc="white", mew=1.4, zorder=3)
    axe.text(1 + i, v + 0.012, f"+{v:.3f}", fontsize=6, ha="center", color=GREEN)
for i, (s, v) in enumerate(amp):
    col = VERM if v > 0.08 else GREY
    mkr = "o" if v > 0.08 else "s"
    axe.plot(4.2 + i, v, mkr, color=col, ms=6, mfc="white", mew=1.4, zorder=3)
    axe.text(4.2 + i, v + 0.012, f"{v:.3f}", fontsize=6, ha="center", color=col)
axe.set_xticks([1, 2, 3, 4.2, 5.2, 6.2])
axe.set_xticklabels(["V1", "V2", "V2_alt", "V1", "V2", "V2_alt*"], fontsize=6)
axe.text(2.0, 0.175, "P-2 midpoint (total effective axis)\nfrozen band ±0.15 dex: 3/3 hit",
         fontsize=6, ha="center", color=GREEN)
axe.text(5.2, 0.175, "P-1b amplitude\nfrozen band 0.08: both arms miss",
         fontsize=6, ha="center", color=VERM)
axe.text(5.2, 0.045, "*ambiguity control,\nnon-scoring (would pass)", fontsize=5.2,
         ha="center", color=GREY)
axe.set_xlim(0.4, 6.8); axe.set_ylim(0.0, 0.20)
axe.set_xlabel("foreground-axis prediction 0.86 µM vs measured 1.24-1.34 µM",
               fontsize=6, color="0.25")
axe.set_ylabel("prediction |error| (dex)")
axe.set_title("Prospective mixed-ligand arm (blind, pre-registered)", fontsize=7)
panel_label(axe, "e")

fig.savefig(OUT + r"\Fig5_receptor_rupture.svg")
fig.savefig(OUT + r"\Fig5_receptor_rupture.png", dpi=200)
print("saved Fig5")
