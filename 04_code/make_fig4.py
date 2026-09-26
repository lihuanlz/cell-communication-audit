# -*- coding: utf-8 -*-
"""Fig. 4 | Blind adjudications bound the claim. Nature-style two-column SVG + 200 dpi PNG.
All numbers transcribed from archived verdict JSONs / SI S7 / Nature_main_v05 legend.
"""
import os
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))

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

C_BLUE = "#0072B2"; C_ORANGE = "#E69F00"; C_GREEN = "#009E73"
C_VERM = "#D55E00"; C_PURPLE = "#CC79A7"; C_SKY = "#56B4E9"
OI7 = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#999999"]

# ---- load P2 units (code 44): 21 WT stratification units, 50min, size, mol_peak ----
csv_path = os.path.join(ROOT, "01_细胞线", "结果", "P2P3P4_盲裁决_本机原件", "代码44_单元表.csv")
rows = list(csv.DictReader(open(csv_path, encoding="utf-8-sig")))
WT = ["ALD3", "DCS2", "DDR2", "HXK1", "RTN2", "SIP18", "TKL2"]
mol21 = [(r["promoter"], r["pair"], float(r["delta"])) for r in rows
         if r["duration"] == "50min" and r["strat"] == "size"
         and r["channel"] == "mol_peak" and r["promoter"] in WT]
assert len(mol21) == 21, len(mol21)
ds = sorted(d for _, _, d in mol21)
print("P2 mol21 range:", ds[0], ds[-1])
time_k3 = {"ALD3": -0.0046, "DCS2": -0.0102, "DDR2": 0.0009, "HXK1": 0.0007,
           "RTN2": -0.0018, "SIP18": -0.0018, "TKL2": 0.0028}  # code44 verdict JSON d_time_k3

fig = plt.figure(figsize=(7.2, 5.2))
gs = fig.add_gridspec(2, 3, width_ratios=[1.55, 1.0, 1.0],
                      left=0.055, right=0.985, top=0.925, bottom=0.105,
                      hspace=0.72, wspace=0.46)

def panel_letter(ax, s, dx=-0.16, dy=1.09):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="top", ha="left")

# ---------------- panel a: P2 Msn2 ----------------
axa = fig.add_subplot(gs[0, 0])
axa.axhspan(-0.01, 0.01, color="0.5", alpha=0.18, zorder=0)
axa.axhline(0.03, color=C_VERM, lw=0.8, ls=(0, (4, 2)), zorder=1)
axa.axhline(0, color="0.4", lw=0.5, zorder=1)
axa.text(20.4, 0.0335, "P2-1 hit line  $\\Delta \\geq +0.03$  (0/7)",
         fontsize=5.4, color=C_VERM, ha="right")
axa.text(4.6, 0.0135, "frozen band $|\\Delta| \\leq 0.01$", fontsize=5.2, color="0.35")
for i, (prom, pair, d) in enumerate(mol21):
    pi = WT.index(prom)
    axa.plot(i, d, "o", ms=3.5, color=OI7[pi], mec="0.25", mew=0.4, zorder=4)
for k in range(1, 7):
    axa.axvline(3 * k - 0.5, color="0.85", lw=0.5, zorder=0)
for k, p in enumerate(WT):
    axa.text(3 * k + 1, -0.068, p, fontsize=4.8, ha="center")
# event-time channel zone (per-promoter d_time, k=3)
x0 = 22.5
axa.axvline(x0 - 1.2, color="0.6", lw=0.6, ls=(0, (2, 2)))
for k, p in enumerate(WT):
    axa.plot(x0 + k, time_k3[p], "s", ms=3.5, color=OI7[WT.index(p)],
             mec="0.25", mew=0.4, zorder=4)
axa.text(x0 + 3.5, 0.092, "event-time $\\Delta_{time}$ ($k=3$)",
         fontsize=5.4, ha="center")
axa.text(x0 + 3.5, 0.081, "6/7 at line;  7/7 ($k$=2), 6/7 ($k$=4)",
         fontsize=5.4, ha="center", color="0.25")
axa.text(x0 + 3.5, -0.068, "same promoter order", fontsize=4.8, ha="center", color="0.4")
axa.text(9.5, 0.068, "molecule-number channel: 21 stratification units",
         fontsize=5.4, ha="center")
axa.text(9.5, 0.057, "(7 promoters x 3 dose pairs),  $\\Delta \\in [-0.049, +0.029]$",
         fontsize=5.4, ha="center", color="0.25")
axa.set_xlim(-0.8, x0 + 6.8); axa.set_ylim(-0.075, 0.100)
axa.set_xticks([])
axa.set_ylabel("stratification penalty $\\Delta$ (AUC)")
axa.set_title("P2: Msn2 (Hansen & Zechner 2021) - amplitude stratification decay",
              fontsize=6.3, pad=16)
panel_letter(axa, "a", dx=-0.09, dy=1.13)

# ---------------- panel b: P3 NF-kB gradient ----------------
axb = fig.add_subplot(gs[0, 1])
groups = [
    ("spatial\n(P3-2)", [(0.907, 0.892, 0.921), (0.819, 0.798, 0.839)], 0.55,
     "hit", C_GREEN),
    ("dose\n(P3-1)", [(0.655, 0.632, 0.678), (0.518, 0.500, 0.538)], 0.60,
     "intermediate", C_ORANGE),
    ("duration\n(P3-3)", [(0.591, 0.546, 0.641), (0.521, 0.481, 0.561)], 0.60,
     "falsified", C_VERM),
]
for gi, (name, pts, line, verdict, vc) in enumerate(groups):
    xc = gi * 2.0
    axb.hlines(line, xc - 0.75, xc + 0.75, color="0.25", lw=0.9, ls=(0, (4, 2)))
    axb.text(xc + 0.80, line - 0.012, f"line {line:.2f}", fontsize=5.0, ha="right",
             va="top", color="0.25")
    for j, (a, lo, hi) in enumerate(pts):
        xx = xc - 0.28 + j * 0.56
        axb.errorbar(xx, a, yerr=[[a - lo], [hi - a]], fmt="o", ms=4, color=vc,
                     mec="0.2", capsize=2, elinewidth=0.8, zorder=5)
        if gi > 0 and j == 1:
            axb.text(xx, lo - 0.020, f"{a:.3f}", fontsize=5.0, ha="center",
                     va="top", color="0.25")
        else:
            axb.text(xx, hi + 0.018, f"{a:.3f}", fontsize=5.0, ha="center",
                     color="0.25")
    axb.text(xc, 0.425, name, fontsize=5.2, ha="center", va="top")
    axb.text(xc, 0.330, verdict, fontsize=5.6, ha="center", color=vc, fontweight="bold")
axb.axhline(0.5, color="0.7", lw=0.5, ls=":")
axb.set_xlim(-1.0, 5.0); axb.set_ylim(0.29, 1.0)
axb.set_xticks([])
axb.set_ylabel("ordinality AUC (95% CI)")
axb.set_title("P3: NF-$\\kappa$B gradient (Son 2022)", fontsize=6.3, pad=16)
panel_letter(axb, "b", dx=-0.30, dy=1.13)

# ---------------- panel c: P4 ERK-KTR ----------------
axc = fig.add_subplot(gs[0, 2])
hist = [0.5819, 0.5633, 0.5636, 0.4695]
uk = [0.5561, 0.5663, 0.5377, 0.5096]
axc.axhline(0.60, color="0.25", lw=0.9, ls=(0, (4, 2)))
axc.text(7.5, 0.605, "frozen line 0.60", fontsize=5.2, ha="right", color="0.25")
axc.axhline(0.50, color="0.6", lw=0.6, ls=":")
axc.text(-0.45, 0.493, "reversal 0.50", fontsize=5.0, ha="left", color="0.45")
xs_h = np.arange(4); xs_u = np.arange(4) + 4.4
axc.plot(xs_h, hist, "o", ms=4, color=C_BLUE, mec="0.2", zorder=5)
axc.plot(xs_u, uk, "s", ms=4, color=C_PURPLE, mec="0.2", zorder=5)
off = [0, -1, 1, -1]
for x, v, o in zip(xs_h, hist, off):
    if o == 0:
        axc.text(x - 0.18, v, f"{v:.4f}", fontsize=4.6, ha="right", va="center",
                 color="0.3")
    else:
        axc.text(x, v + 0.013 * o + (0.004 if o > 0 else -0.006), f"{v:.4f}",
                 fontsize=4.6, ha="center", color="0.3")
for x, v, o in zip(xs_u, uk, off):
    if o == 0:
        axc.text(x - 0.32, v + 0.018, f"{v:.4f}", fontsize=4.6, ha="right",
                 va="center", color="0.3", zorder=6)
    else:
        axc.text(x, v + 0.013 * o + (0.004 if o > 0 else -0.006), f"{v:.4f}",
                 fontsize=4.6, ha="center", color="0.3", zorder=6)
axc.text(1.5, 0.445, "histamine\n$n$ = 7,393 cells", fontsize=5.2, ha="center",
         va="top", color=C_BLUE)
axc.text(5.9, 0.445, "UK14304\n$n$ = 7,508 cells", fontsize=5.2, ha="center",
         va="top", color=C_PURPLE)
axc.text(3.7, 0.695, "all 8 AUCs below line (max 0.5819): falsified\n"
         "stratification: 5/10 units $|\\Delta| \\geq 0.03$ (falsified)",
         fontsize=5.2, ha="center", va="top")
axc.set_xlim(-0.6, 7.9); axc.set_ylim(0.40, 0.70)
axc.set_xticks([])
axc.set_ylabel("event-time decoding AUC")
axc.set_title("P4: ERK-KTR (Chavez-Abiega 2022)", fontsize=6.3, pad=16)
panel_letter(axc, "c", dx=-0.30, dy=1.13)

# ---------------- panel d: P5 GPCR -> Ca2+ ----------------
axd = fig.add_subplot(gs[1, 0])
axd.add_patch(Rectangle((-0.7, 0.40), 4.5, 0.20, color="0.5", alpha=0.18, zorder=0))
tau_auc = [("D1-D2", 0.629), ("D2-D3", 0.659), ("D5-D6", 0.568), ("D6-D7", 0.520)]
xs = [0, 1.25, 2.5, 3.75]
for x, (lab, v) in zip(xs, tau_auc):
    out = not (0.40 <= v <= 0.60)
    axd.plot(x, v, "o", ms=4.5, color=C_VERM if out else C_BLUE, mec="0.2", zorder=5)
    axd.text(x, v + 0.025, f"{v:.3f}", fontsize=5.0, ha="center", color="0.3")
    axd.text(x, 0.055, lab, fontsize=4.8, ha="center")
axd.text(1.875, 0.86, "dead-zone clause: band [0.40, 0.60]\n"
         "2/4 units out of band: falsified", fontsize=5.4, ha="center", va="top")
axd.text(1.875, 0.455, "frozen band", fontsize=4.8, ha="center", color="0.4")
# right zone: population counting arm (descriptive)
rate = [0.015, 0.272, 0.544, 0.769, 0.826, 0.892, 0.918]
xr = np.arange(7) * 1.05 + 5.6
axd.axvline(4.7, color="0.6", lw=0.6, ls=(0, (2, 2)))
axd.plot(xr, rate, "-o", ms=3.6, color=C_GREEN, mec="0.2", lw=1.1, zorder=5)
for x, v in zip(xr, rate):
    if v < 0.05:
        axd.text(x - 0.12, v + 0.01, f"{v:.3f}", fontsize=4.6, ha="right",
                 va="center", color="0.3")
    else:
        axd.text(x, v + 0.032, f"{v:.3f}", fontsize=4.6, ha="center", color="0.3")
axd.text(8.75, 0.26, "population counting arm (descriptive):\n"
         "response rate 0.015 -> 0.918,\nmonotonic over D1-D7",
         fontsize=5.4, ha="center", va="top")
for x, lab in zip(xr, ["D1", "D2", "D3", "D4", "D5", "D6", "D7"]):
    axd.text(x, 0.055, lab, fontsize=4.8, ha="center")
axd.set_xlim(-0.9, 12.7); axd.set_ylim(0.0, 1.02)
axd.set_xticks([])
axd.set_ylabel("AUC / response rate")
axd.set_title("P5: GPCR-to-Ca$^{2+}$ boundary predictor (Keshelava 2018)",
              fontsize=6.3, pad=4)
panel_letter(axd, "d", dx=-0.09, dy=1.09)

# ---------------- panel e: P6/P7 forest plot ----------------
axe = fig.add_subplot(gs[1, 1:])
forest = [
    ("P6  NF-$\\kappa$B sequential\nstim. (Wang 2022)", -0.0014, -0.2518, 0.1533, C_BLUE),
    ("P7  chemotaxis ladder\n(Moore 2024)", 0.1181, -0.2956, 0.4677, C_GREEN),
    ("P7  erratum arm\n(49b re-adjudication)", -0.0909, -0.4446, 0.2766, C_ORANGE),
]
axe.axvspan(-0.15, 0.15, color="0.5", alpha=0.15, zorder=0)
for ln in (-0.15, 0.15):
    axe.axvline(ln, color="0.25", lw=0.8, ls=(0, (4, 2)))
axe.axvline(0, color="0.4", lw=0.6)
axe.text(0.0, 3.30, "decision lines $\\pm 0.15$", fontsize=5.4, color="0.25",
         ha="center")
for i, (lab, d, lo, hi, cc) in enumerate(forest):
    y = 2 - i
    axe.plot([lo, hi], [y, y], color=cc, lw=1.4, solid_capstyle="round")
    axe.plot([lo, lo], [y - 0.09, y + 0.09], color=cc, lw=1.0)
    axe.plot([hi, hi], [y - 0.09, y + 0.09], color=cc, lw=1.0)
    axe.plot(d, y, "o", ms=5, color=cc, mec="0.2", zorder=5)
    axe.text(-0.705, y, lab, fontsize=5.2, ha="left", va="center")
    axe.text(0.0, y - 0.34, f"$\\Delta R^2$ = {d:+.4f}  [{lo:+.3f}, {hi:+.3f}]",
             fontsize=5.2, ha="center", va="center", color="0.25")
axe.text(0.0, -0.72, "all CIs cross zero: power-limited intermediates",
         fontsize=5.8, ha="center", color="0.25")
axe.set_xlim(-0.72, 0.60); axe.set_ylim(-0.95, 3.55)
axe.set_yticks([])
axe.set_xticks([-0.4, -0.2, 0, 0.2, 0.4])
axe.set_xlabel(r"$\Delta R^2$ (fold-change vs absolute-dose attribution)")
axe.set_title("P6/P7: fold-change attribution", fontsize=6.3, pad=4)
panel_letter(axe, "e", dx=-0.08, dy=1.09)

fig.savefig(os.path.join(HERE, "Fig4_adjudications.svg"))
fig.savefig(os.path.join(HERE, "Fig4_adjudications.png"), dpi=200)
print("Fig4 written")
