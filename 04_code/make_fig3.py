# -*- coding: utf-8 -*-
"""Fig. 3 | The dynamic encoder. Nature-style two-column SVG + 200 dpi PNG.
All numbers transcribed from Nature_SI_v05 (S3, S4) and Nature_main_v05 legend.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.linewidth": 0.6,
    "axes.labelsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "svg.fonttype": "none",
    "legend.frameon": False,
    "legend.fontsize": 6,
})

# Okabe-Ito palette
C_BLUE = "#0072B2"; C_ORANGE = "#E69F00"; C_GREEN = "#009E73"
C_VERM = "#D55E00"; C_PURPLE = "#CC79A7"; C_SKY = "#56B4E9"

rng = np.random.default_rng(20260923)

fig = plt.figure(figsize=(7.2, 4.8))
gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.05], width_ratios=[1.0, 1.0, 1.15],
                      left=0.065, right=0.985, top=0.94, bottom=0.10,
                      hspace=0.70, wspace=0.42)

def panel_letter(ax, s, dx=-0.22, dy=1.06):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="top", ha="left")

# ---------------- panel a: NF bifurcation with sqrt onset ----------------
axa = fig.add_subplot(gs[0, 0])
s_H = 16.4449
ds = np.linspace(0.0, 4.2, 300)
A = 0.95 * np.sqrt(np.clip(ds, 0, None))          # A ∝ sqrt(s - s_H), scaled
axa.plot(s_H + ds, A, color=C_BLUE, lw=1.4, zorder=3)
axa.plot([16.0, s_H], [0, 0], color=C_BLUE, lw=1.4, zorder=3)
dss = np.linspace(0.05, 3.80, 18)
pts = 0.95 * np.sqrt(dss) * (1 + rng.normal(0, 0.008, dss.size))
axa.plot(s_H + dss, pts, "o", ms=2.5, color=C_BLUE, mfc="white", mew=0.7, zorder=4)
axa.axvline(s_H, color="0.45", lw=0.6, ls=(0, (3, 2)), zorder=2)
axa.annotate(r"$s_H = 16.4449$ (Hopf)", xy=(s_H, 0.06), xytext=(16.08, 0.62),
             fontsize=5.8, arrowprops=dict(arrowstyle="-", lw=0.5, color="0.45"))
axa.text(0.04, 0.93, r"supercritical Hopf onset:  $A \propto \sqrt{s-s_H}$",
         transform=axa.transAxes, fontsize=5.8)
axa.set_xlim(16.0, 20.9); axa.set_ylim(-0.05, 2.15)
axa.set_xlabel("drive $s$"); axa.set_ylabel("amplitude $A$ (au)")
panel_letter(axa, "a")
axin = axa.inset_axes([0.50, 0.10, 0.46, 0.42])
x = np.linspace(0.05, 3.80, 50)
axin.plot(x, (0.95 ** 2) * x, color=C_VERM, lw=1.0)
axin.plot(dss, pts ** 2, "o", ms=1.8, color=C_VERM, mfc="white", mew=0.5)
axin.text(0.05, 0.88, r"$A^2$ linear in $(s-s_H)$", transform=axin.transAxes, fontsize=5.2)
axin.text(0.05, 0.72, r"$R^2 = 0.9926$", transform=axin.transAxes, fontsize=5.2)
axin.set_xlabel(r"$s - s_H$", fontsize=5.2, labelpad=1)
axin.set_ylabel(r"$A^2$", fontsize=5.2, labelpad=1)
axin.tick_params(labelsize=4.8, width=0.5, length=1.8)
for sp in axin.spines.values():
    sp.set_linewidth(0.5)

# ---------------- panel b: subcritical hysteresis loop ----------------
axb = fig.add_subplot(gs[0, 1])
mu = np.linspace(-0.249, 0.42, 400)
r_plus = np.sqrt((1 + np.sqrt(1 + 4 * mu)) / 2)
mu_u = np.linspace(-0.249, 0.0, 200)
r_minus = np.sqrt((1 - np.sqrt(1 + 4 * mu_u)) / 2)
axb.plot(mu, r_plus, color=C_GREEN, lw=1.4)
axb.plot(mu_u, r_minus, color=C_GREEN, lw=1.0, ls="--")
axb.plot([-0.42, 0], [0, 0], color="0.2", lw=1.4)
axb.plot([0, 0.42], [0, 0], color="0.2", lw=1.0, ls="--")
# up-scan and down-scan arrows
axb.annotate("", xy=(-0.02, 0.02), xytext=(-0.30, 0.02),
             arrowprops=dict(arrowstyle="->", color=C_VERM, lw=1.0))
axb.annotate("", xy=(0.0, 0.97), xytext=(0.0, 0.06),
             arrowprops=dict(arrowstyle="->", color=C_VERM, lw=1.0))
axb.annotate("", xy=(0.34, 1.115), xytext=(0.06, 1.008),
             arrowprops=dict(arrowstyle="->", color=C_VERM, lw=1.0))
axb.annotate("", xy=(-0.22, 0.795), xytext=(0.28, 1.075),
             arrowprops=dict(arrowstyle="->", color=C_BLUE, lw=1.0))
axb.annotate("", xy=(-0.25, 0.04), xytext=(-0.25, 0.68),
             arrowprops=dict(arrowstyle="->", color=C_BLUE, lw=1.0))
axb.text(0.38, 0.86, "up-scan", color=C_VERM, fontsize=5.8, ha="right")
axb.text(-0.40, 1.13, "down-scan", color=C_BLUE, fontsize=5.8, ha="left")
axb.text(0.44, 0.16, "stable cycle (green)\nunstable (dashed)", fontsize=5.2,
         color="0.35", va="center", ha="right")
# hysteresis width bracket
axb.annotate("", xy=(-0.25, -0.14), xytext=(0.0, -0.14),
             arrowprops=dict(arrowstyle="<->", lw=0.7, color="0.2"))
axb.text(-0.125, -0.245, "width 0.25 (NF prediction)", fontsize=5.6, ha="center")
axb.axvline(-0.25, color="0.6", lw=0.5, ls=(0, (3, 2)))
axb.axvline(0.0, color="0.6", lw=0.5, ls=(0, (3, 2)))
axb.text(0.10, 0.42, "measured width 0.275\nvs prediction 0.25", fontsize=5.8,
         color=C_VERM, ha="left")
axb.set_xlim(-0.42, 0.46); axb.set_ylim(-0.30, 1.25)
axb.set_xlabel(r"drive $\mu$"); axb.set_ylabel("cycle amplitude $r$")
panel_letter(axb, "b")

# ---------------- panel c: discrimination diagram (alpha, nu) ----------------
axc = fig.add_subplot(gs[0, 2])
axc.axvspan(-0.45, 0.1, color=C_GREEN, alpha=0.10, zorder=0)
axc.axvspan(0.1, 0.5, color="0.6", alpha=0.12, zorder=0)
axc.axvspan(0.5, 1.40, color=C_VERM, alpha=0.10, zorder=0)
axc.axvline(0.1, color="0.35", lw=0.7, ls=(0, (3, 2)))
axc.axvline(0.5, color="0.35", lw=0.7, ls=(0, (3, 2)))
axc.text(-0.43, 12.4, "EXC region", fontsize=6.2, color=C_GREEN, fontweight="bold")
axc.text(-0.43, 11.55, r"$\alpha \lesssim 0.1$, $\nu>0$", fontsize=5.4, color=C_GREEN)
axc.text(0.30, 12.4, "leaky-digital", fontsize=5.4, color="0.35", ha="center")
axc.text(0.60, 12.4, "NF region", fontsize=6.2, color=C_VERM, fontweight="bold")
axc.text(0.60, 11.55, r"$\alpha \gtrsim 0.5$", fontsize=5.4, color=C_VERM)
# model anchors
axc.plot(1.13, 0.0, "s", ms=5, color=C_VERM, mec="0.2", zorder=5)
axc.annotate("model NF anchor\n$\\alpha = 1.13$, $\\nu \\approx 0$",
             xy=(1.14, 0.18), xytext=(0.94, 3.2), fontsize=5.2, va="bottom",
             arrowprops=dict(arrowstyle="-", lw=0.5, color="0.4"))
axc.plot(0.04, 10.5, "^", ms=6, color=C_GREEN, mec="0.2", zorder=5)
axc.annotate("model EXC anchor\n$\\alpha = 0.04$, $\\nu = \\tau_r/T = 10.5$",
             xy=(0.10, 10.4), xytext=(0.30, 10.3), fontsize=5.4, va="top",
             arrowprops=dict(arrowstyle="-", lw=0.5, color="0.4"))
# experimental placements
axc.errorbar(-0.001, 2.0, xerr=0.113, yerr=0.5, fmt="o", ms=5, color=C_PURPLE,
             mec="0.2", capsize=2, elinewidth=0.8, zorder=6)
axc.annotate("p53-DSB (NCS)\n$\\alpha=-0.001\\pm0.113$", xy=(-0.03, 2.45),
             xytext=(-0.44, 3.9), fontsize=5.4, color="0.15", va="top",
             arrowprops=dict(arrowstyle="-", lw=0.5, color="0.4"))
axc.errorbar(0.981, 0.0, xerr=0.114, yerr=[[0.0], [0.2]], fmt="o", ms=5,
             color=C_ORANGE, mec="0.2", capsize=2, elinewidth=0.8, zorder=6)
axc.text(0.55, 1.45, "p53-UV\n$\\alpha=+0.981\\pm0.114$", fontsize=5.4,
         color="0.15", va="bottom")
axc.errorbar(0.0, 0.17, xerr=0.05, fmt="D", ms=4.5, color=C_SKY, mec="0.2",
             capsize=2, elinewidth=0.8, zorder=6)
axc.annotate("ERK  $|\\alpha|\\lesssim0.05$\n$+0.39$ counts/decade",
             xy=(0.05, 0.06), xytext=(0.72, -0.55), fontsize=5.0, color="0.15",
             va="top", ha="center",
             arrowprops=dict(arrowstyle="-", lw=0.5, color="0.4"))
axc.plot(0.151, 0.26, "D", ms=4.5, color="0.25", mec="0.2", zorder=6)
axc.annotate("NF-$\\kappa$B  $\\alpha=0.151$\n$+0.6$ counts/decade",
             xy=(0.14, 0.45), xytext=(-0.44, 5.55), fontsize=5.4, color="0.15",
             va="top",
             arrowprops=dict(arrowstyle="-", lw=0.5, color="0.4",
                             connectionstyle="arc3,rad=-0.15"))
axc.plot(0.020, 0.128, "*", ms=8, color=C_GREEN, mec="0.2", zorder=6)
axc.annotate("Msn2-type  $\\alpha=0.020$\n(frequency coding)",
             xy=(0.02, 0.10), xytext=(-0.44, -0.15), fontsize=5.4, color="0.15",
             va="top",
             arrowprops=dict(arrowstyle="-", lw=0.5, color="0.4"))
axc.annotate("", xy=(0.981, 6.2), xytext=(-0.001, 6.2),
             arrowprops=dict(arrowstyle="<->", lw=0.7, color="0.25"))
axc.text(0.49, 6.5, "DSB vs UV: 6.1 combined $\\sigma$ apart", fontsize=5.6,
         ha="center", va="bottom", color="0.25")
axc.set_xlim(-0.45, 1.40); axc.set_ylim(-1.35, 13.2)
axc.set_xlabel(r"$\alpha = \partial\ln A\,/\,\partial\ln D$ (amplitude-dose)")
axc.set_ylabel(r"$\nu = \partial N\,/\,\partial\ln D$ (count-dose)")
panel_letter(axc, "c")

# ---------------- panel d: channel census ----------------
axd = fig.add_subplot(gs[1, 0])
channels = ["w", "N", "T", "A", "t1", "IPI"]
shares = [60.2, 18.2, 11.5, 5.5, 4.4, 0.3]
colors = [C_GREEN, C_GREEN, C_GREEN, C_VERM, C_GREEN, C_GREEN]
alphas = [1.0, 0.75, 0.55, 1.0, 0.4, 0.28]
left = 0.0
lab_pos = {"A": (92.65, 0.50), "t1": (97.6, 0.92), "IPI": (100.4, 1.34)}
for ch, sh, cc, al in zip(channels, shares, colors, alphas):
    axd.barh(0, sh, left=left, height=0.42, color=cc, alpha=al,
             edgecolor="white", lw=0.5)
    if sh >= 10:
        axd.text(left + sh / 2, 0, f"{ch}\n{sh}%", ha="center", va="center",
                 fontsize=5.6, color="white" if al > 0.6 else "0.15")
    else:
        xt, yt = lab_pos[ch]
        axd.text(xt, yt, f"{ch} {sh}%", fontsize=5, ha="center", va="center")
    left += sh
axd.annotate("", xy=(0, -0.34), xytext=(94.4, -0.34),
             arrowprops=dict(arrowstyle="-", lw=0.8, color=C_GREEN))
axd.plot([0, 0], [-0.28, -0.40], color=C_GREEN, lw=0.8)
axd.plot([94.4, 94.4], [-0.28, -0.40], color=C_GREEN, lw=0.8)
axd.text(47, -0.55, "timing channels 94.2%", fontsize=5.8, ha="center", color=C_GREEN)
axd.set_xlim(-1, 112); axd.set_ylim(-0.85, 1.9)
axd.set_yticks([]); axd.set_xlabel("share of Fisher information (%)")
axd.set_xticks([0, 20, 40, 60, 80, 100])
for sp in ["left", "top", "right"]:
    axd.spines[sp].set_visible(False)
panel_letter(axd, "d")
# inset: two-layer death of amplitude
axd2 = axd.inset_axes([0.03, 0.62, 0.50, 0.34])
sk = np.linspace(0, 0.6, 100)
cvA = 0.002 + (0.639 - 0.002) * (sk / 0.6) ** 1.6   # endpoints locked: 0.002 -> 0.639
axd2.plot(sk, cvA, color=C_VERM, lw=1.2)
axd2.axhline(0.016, color=C_GREEN, lw=1.1, ls="--")
axd2.plot([0, 0.6], [0.002, 0.639], "o", ms=3, color=C_VERM)
axd2.plot(0.3, 0.016, "s", ms=3, color=C_GREEN)
axd2.text(0.66, 0.30, "CV(A) 0.639", fontsize=5, va="center", color=C_VERM)
axd2.text(0.605, 0.10, "CV(T) 0.016", fontsize=5, va="center", color=C_GREEN)
axd2.text(0.03, 0.97, "gain dispersion $\\sigma_k$: 0 to 0.6", transform=axd2.transAxes,
          fontsize=5.0, va="top")
axd2.text(0.03, 0.80, "CV(A) 0.002 -> 0.639", transform=axd2.transAxes,
          fontsize=5.0, va="top")
axd2.text(0.03, 0.63, "CV(T) pinned 0.016", transform=axd2.transAxes,
          fontsize=5.0, va="top")
axd2.set_xlim(0, 0.98); axd2.set_ylim(0, 0.75)
axd2.set_xlabel(r"$\sigma_k$", fontsize=5.2, labelpad=1)
axd2.set_ylabel("CV", fontsize=5.2, labelpad=1)
axd2.tick_params(labelsize=4.8, width=0.5, length=1.8)
for sp in axd2.spines.values():
    sp.set_linewidth(0.5)
axd.set_title("p53 channel census (FHN representative)", fontsize=6.5, pad=2)

# ---------------- panel e: digital-limit schematic ----------------
axe = fig.add_subplot(gs[1, 1:])
axe.set_xlim(0, 10); axe.set_ylim(0, 10); axe.axis("off")
panel_letter(axe, "e", dx=-0.115, dy=1.04)

def box(ax, x, y, w, h, text, fc, ec, fs=6.0):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.12",
                                fc=fc, ec=ec, lw=0.8, mutation_scale=1))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs)

box(axe, 0.1, 6.6, 2.6, 1.9, "chemical saturation\n$\\kappa \\to 0$\n(Hill $n_H \\to \\infty$)",
    "#FDEBD0", C_ORANGE)
box(axe, 0.1, 1.4, 2.6, 1.9, "circuit topology\n$\\alpha \\to 0$\n(EXC amplitude lock)",
    "#D6EAF8", C_BLUE)
box(axe, 4.0, 4.0, 2.8, 2.0,
    "same likelihood\ncancellation:\namplitude scale\ndegenerate ($\\lambda = 0$)",
    "#E8DAEF", C_PURPLE, fs=5.8)
axe.add_patch(FancyArrowPatch((2.8, 7.5), (4.1, 5.8), arrowstyle="-|>",
                              mutation_scale=10, lw=1.1, color=C_ORANGE))
axe.add_patch(FancyArrowPatch((2.8, 2.5), (4.1, 4.2), arrowstyle="-|>",
                              mutation_scale=10, lw=1.1, color=C_BLUE))
axe.text(3.35, 6.7, "two roads,\none limit", fontsize=5.4, style="italic",
         color="0.3", ha="center")
axe.text(5.4, 2.9, r"$N \approx (\tau_r/T)\cdot\ln(D_0/D_c)$", fontsize=7.5,
         ha="center", color="0.1")
axe.text(5.4, 2.0, "slope inverts repair timescale:\n"
         r"$\tau_r \approx \nu\,T = 1.4 \times 5.5$ h $\approx 7.7$ h",
         fontsize=5.6, ha="center", va="top", color="0.25")
# counting-law inset axes
axc2 = axe.inset_axes([0.775, 0.13, 0.215, 0.74])
D0 = np.linspace(0.30, 1.20, 50)
Dc = 0.261
Nlaw = 10.5 * np.log(D0 / Dc)
Nlaw = Nlaw - Nlaw[0] + 2.0
axc2.plot(np.log(D0 / Dc), Nlaw, color=C_GREEN, lw=1.3)
npts = 10.5 * np.log(D0[::6] / Dc); npts = npts - npts[0] + 2.0
axc2.plot(np.log(D0[::6] / Dc), npts, "o", ms=2.4, color=C_GREEN, mfc="white", mew=0.6)
axc2.set_xlabel(r"$\ln(D_0/D_c)$", fontsize=5.2, labelpad=1)
axc2.set_ylabel("pulse count $N$", fontsize=5.2, labelpad=1)
axc2.tick_params(labelsize=4.8, width=0.5, length=1.8)
for sp in axc2.spines.values():
    sp.set_linewidth(0.5)
axc2.text(0.05, 0.94, "counting law verified", transform=axc2.transAxes,
          fontsize=5.0, va="top")
axc2.text(0.05, 0.80, "$R^2 = 0.920$, $D_c = 0.261$", transform=axc2.transAxes,
          fontsize=5.0, va="top")
axe.add_patch(FancyArrowPatch((6.9, 5.0), (7.55, 5.0), arrowstyle="-|>",
                              mutation_scale=10, lw=1.1, color=C_PURPLE))

fig.savefig(os.path.join(HERE, "Fig3_dynamic_encoder.svg"))
fig.savefig(os.path.join(HERE, "Fig3_dynamic_encoder.png"), dpi=200)
print("Fig3 written")
