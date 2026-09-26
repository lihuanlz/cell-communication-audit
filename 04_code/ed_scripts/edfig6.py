import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from edstyle import setup, panel_label, schematic_tag, save, OI
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

setup()
ROOT = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911"
OUT = os.path.join(ROOT, r"05_主线纲领与设计\论文_细胞通讯审计_2026-09-23\figures_svg\ED")

fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.9))
axs = axs.ravel()

# ---- a: degeneracy hierarchy G0 < G1 < G2 ----
ax = axs[0]
ax.axis("off")
boxes = [
    (0.03, 0.06, 0.94, 0.66, "G2: monotone reparameterization", OI["verm"]),
    (0.10, 0.14, 0.80, 0.48, "G1: additive baseline", OI["orange"]),
    (0.17, 0.22, 0.66, 0.30, "G0: multiplicative", OI["blue"]),
]
for x, y, w, h, lab, c in boxes:
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.008",
                                fc=c, ec=c, alpha=0.14, lw=0.9))
    ax.text(x + 0.03, y + h - 0.055, lab, fontsize=6, color=c, fontweight="bold", va="top")
ax.text(0.20, 0.30, "counting statistics survive", fontsize=5.4, ha="left", color="0.25")
ax.text(0.13, 0.185, "ratio statistics die", fontsize=5.4, color="0.25")
ax.text(0.06, 0.095, "event time / order survive;  L1' max invariant R_y", fontsize=5.4, color="0.25")
ax.text(0.5, 0.82, "nested degeneracy groups: G0 < G1 < G2", fontsize=6.2, ha="center", color="0.15")
ax.set_xlim(0, 1); ax.set_ylim(0, 0.9)
ax.set_title("digital-limit hierarchy: what survives each group")
schematic_tag(ax)

# ---- b: code13 verification spectrum ----
ax = axs[1]
xg = np.arange(4)
w = 0.36
pure = [26935, 1898, 0.0, 0.0]
withamp = [106667, 26935, 1898, 0.0]
ax.bar(xg - w / 2, [v if v > 0 else 0.9 for v in pure], w, color=OI["blue"], label="counting only")
ax.bar(xg + w / 2, [v if v > 0 else 0.9 for v in withamp], w, color=OI["verm"], label="counting + amplitude")
for xi, v in zip(xg - w / 2, pure):
    ax.text(xi, v * 1.15 if v > 0 else 1.3, f"{v:g}" if v > 0 else "0 (exact)", ha="center", fontsize=5.2, color=OI["blue"])
for xi, v in zip(xg + w / 2, withamp):
    ax.text(xi, v * 1.15 if v > 0 else 1.3, f"{v:g}" if v > 0 else "0 (exact)", ha="center", fontsize=5.2, color=OI["verm"])
ax.set_yscale("log"); ax.set_ylim(0.8, 3e6)
ax.set_xticks(xg); ax.set_xticklabels([f"lam_{i+1}" for i in xg], fontsize=6)
ax.set_ylabel("Fisher eigenvalue")
ax.text(0.02, 0.04, "zero direction (0,0,1,-1)/sqrt(2): cosine 1.000000\n"
        "CRB ratios 1.24 / 1.36 / 0.95;  structure block kappa = 14.2",
        transform=ax.transAxes, fontsize=5.2, va="bottom", color="0.25")
ax.set_title("code 13: exact zero direction survives amplitude channel")
ax.legend(fontsize=5.2, loc="upper right")
schematic_tag(ax)

# ---- c: noise concentration bound N1 + numerical pins ----
ax = axs[2]
t = np.linspace(0, 0.6, 200)
for n in [50, 200, 1000]:
    ax.plot(t, 2 * np.exp(-n * t ** 2 / 2), "-", lw=1.0, label=f"n = {n}")
ax.set_yscale("log")
ax.set_xlabel("deviation t"); ax.set_ylabel("tail bound 2 exp(-n t$^2$/2)")
ax.set_title("noise concentration (N1) and numerical pins")
ax.legend(fontsize=5.2, loc="upper right")
ax.text(0.03, 0.30, "theory 0.053671 vs empirical 0.053472\n"
        "40% plateau: total distortion 0.1840 vs 0.1843\n"
        "p53 anchor: Phi(-4.71) ~ 1.2e-6",
        transform=ax.transAxes, fontsize=5.2, color="0.25")
schematic_tag(ax)

# ---- d: conjecture decision tree + code25 verdict ----
ax = axs[3]
ax.axis("off")
def box(x, y, w, h, text, fc, tc="0.1", fs=5.4):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.01", fc=fc, ec="0.4", lw=0.7, alpha=0.9))
    ax.text(x + w / 2, y + h / 2, text, fontsize=fs, ha="center", va="center", color=tc)
def arr(x1, y1, x2, y2, lab=None, dx=0.02):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="->", lw=0.8, color="0.35"))
    if lab:
        ax.text((x1 + x2) / 2 + dx, (y1 + y2) / 2, lab, fontsize=5, color="0.4")
box(0.30, 0.80, 0.40, 0.13, "cellular statistic\nunder audit", "#DCE9F5")
box(0.05, 0.50, 0.42, 0.16, "event time / order?\n(invariant under G2)", "#E4F0E4")
box(0.55, 0.50, 0.42, 0.16, "amplitude / ratio?\n(dies under G1, G2)", "#F6E3DE")
arr(0.42, 0.80, 0.26, 0.66)
arr(0.58, 0.80, 0.76, 0.66)
box(0.05, 0.20, 0.42, 0.16, "admissible currency\n(L1': R_y max invariant)", "#CFE6CF")
box(0.55, 0.20, 0.42, 0.16, "not identifiable:\ndrop or fix scale", "#F0CFC7")
arr(0.26, 0.50, 0.26, 0.36)
arr(0.76, 0.50, 0.76, 0.36)
ax.text(0.5, 0.06, "code 25 first test (Morris-Lecar): 2/3 passed -\n"
        "cosine 0.0647 fails; conjecture remains OPEN",
        fontsize=5.6, ha="center", color=OI["verm"])
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_title("admissibility decision tree (conjecture status)")
schematic_tag(ax)

# ---- e: energy cost crossover ----
ax = axs[4]
L = np.logspace(0, 1.6, 200)
alpha = 7.0
EA = 2 * L ** 2   # 2 kT L^2 D^2 (units kT D^2)
EB = alpha * L    # alpha kT L D^2
ax.plot(L, EA, "-", color=OI["blue"], lw=1.2, label="E_A = 2 kT L$^2$ D$^2$ (counting)")
ax.plot(L, EB, "-", color=OI["verm"], lw=1.2, label="E_B = alpha kT L D$^2$ (amplitude)")
ax.axvline(alpha / 2, color="0.4", ls=":", lw=0.7)
ax.text(3.65, 60, "L* = alpha/2 ~ 3-4\ncorrelation error sqrt(alpha) ~ 2.6", fontsize=5.2, color="0.3")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("dynamic range L"); ax.set_ylabel("energy cost (kT D$^2$)")
ax.set_title("energy crossover favors counting at large L")
ax.legend(fontsize=5.2, loc="upper left")
schematic_tag(ax)

# ---- f: information ceilings ----
ax = axs[5]
names = ["simulated channel\n(eps0 = 0.2)", "counting soft ceiling\n(D0/Dc = 100)"]
vals = [2.35, 4.5]
bars = ax.barh(names, vals, color=[OI["sky"], OI["blue"]], height=0.5)
for b, v in zip(bars, vals):
    ax.text(v + 0.08, b.get_y() + b.get_height() / 2, f"{v:.2f} bits", va="center", fontsize=6)
ax.set_xlim(0, 5.6)
ax.set_xlabel("information ceiling (bits)")
ax.set_title("channel capacity ceilings (archived)")
ax.tick_params(axis="y", labelsize=6)
schematic_tag(ax)

for i, ax in enumerate(axs):
    panel_label(ax, "abcdef"[i])

fig.suptitle("EDFig. 6  Digital-limit theorems: degeneracy hierarchy, Fisher zero direction, noise concentration,\n"
             "admissibility decision tree, energy crossover, information ceilings (schematic redraws of archived results)",
             fontsize=8)
fig.tight_layout(rect=[0, 0, 1, 0.915])
save(fig, os.path.join(OUT, "EDFig6_digital_limit_theorems"))
