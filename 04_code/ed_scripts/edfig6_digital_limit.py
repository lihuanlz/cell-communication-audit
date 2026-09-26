"""EDFig. 6: Digital-limit theorem diagrams and conjecture decision tree.

Panels:
  a  two-roads theorem diagram (S5.1) + S5.6 energy crossover inset
  b  Theorem T' stratified classification G0 < G1 < G2 (S5.2) + numeric pins
  c  conjecture decision tree (S1.5 legs, S5.4 conjecture, Morris-Lecar open)
  d  Theorem N noise-robustness sketch (S5.2, eqs S5.1/S5.2)
All numbers transcribed from _archive_旧版本/Nature_SI_v04_2026-09-23.md, section S5.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from edstyle import setup, panel_label, schematic_tag, save, OI
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

setup()
ROOT = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911"
OUT = os.path.join(ROOT, r"05_主线纲领与设计\论文_细胞通讯审计_2026-09-23\figures_svg\ED")

GREY = "0.45"

fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.4))

def box(ax, x, y, w, h, text, fc, ec=None, tc="0.1", fs=5.4, lw=0.7, bold=False):
    ec = ec if ec is not None else "0.4"
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.008",
                                fc=fc, ec=ec, lw=lw))
    ax.text(x + w / 2, y + h / 2, text, fontsize=fs, ha="center", va="center",
            color=tc, fontweight="bold" if bold else "normal")

def arr(ax, x1, y1, x2, y2, color="0.35"):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="->", lw=0.8, color=color))

# ================= panel a: two-roads theorem diagram =================
ax = axs[0, 0]
ax.axis("off")
ax.set_xlim(0, 1); ax.set_ylim(0, 1)

ax.text(0.5, 0.985, "two constructive roads to the counting likelihood (S5.1)",
        fontsize=6.4, ha="center", va="top", color="0.15", fontweight="bold")

box(ax, 0.02, 0.72, 0.30, 0.20,
    "PdPC static sensor\nchemical saturation\nkappa -> 0",
    "#D6E6F2", ec=OI["blue"], tc=OI["blue"])
box(ax, 0.68, 0.72, 0.30, 0.20,
    "p53 pulse train\ncircuit topology\nalpha -> 0",
    "#F6E0D6", ec=OI["verm"], tc=OI["verm"])

box(ax, 0.02, 0.46, 0.30, 0.16,
    "independent Bernoulli\nmolecules, hence binomial",
    "#EAF1F8", ec=OI["blue"], fs=5.0)
box(ax, 0.68, 0.46, 0.30, 0.16,
    "Poisson count jitter\nafter scale cancellation",
    "#FAEDE6", ec=OI["verm"], fs=5.0)

arr(ax, 0.17, 0.72, 0.17, 0.62, color=OI["blue"])
arr(ax, 0.83, 0.72, 0.83, 0.62, color=OI["verm"])
arr(ax, 0.17, 0.46, 0.42, 0.33, color=OI["blue"])
arr(ax, 0.83, 0.46, 0.58, 0.33, color=OI["verm"])

box(ax, 0.26, 0.12, 0.40, 0.22,
    "shared counting likelihood\n"
    "p(n) = [N! / prod n_k!] prod q_k^n_k\n"
    "n sufficient (C2); invariant +\n"
    "relevant => counting type (C3, S2.4)",
    "#DDEBDD", ec=OI["green"], fs=4.9, bold=True)

ax.text(0.5, 0.665, "conditional theorem: H1-H3 in, C1-C3 out\n"
        "constructive proof on two circuits; CV excluded by tightened C3",
        fontsize=4.8, ha="center", va="center", color="0.3")

# inset: S5.6 energy crossover
axi = ax.inset_axes([0.015, 0.045, 0.235, 0.27])
L = np.logspace(0, 1.3, 200)
alpha = 7.0
axi.plot(L, 2 * L ** 2, "-", color=OI["blue"], lw=1.0,
         label="E_A = 2kTL$^2$D$^2$")
axi.plot(L, alpha * L, "-", color=OI["green"], lw=1.0,
         label="E_B = alpha kTLD$^2$")
axi.axvline(alpha / 2, color="0.4", ls=":", lw=0.7)
axi.set_xscale("log"); axi.set_yscale("log")
axi.set_xlim(1, 20); axi.set_ylim(1.5, 1e3)
axi.set_xticks([1, 10]); axi.set_yticks([1e1, 1e3])
axi.tick_params(labelsize=4.2, length=1.5, width=0.5)
for s in axi.spines.values():
    s.set_linewidth(0.5)
axi.set_title("S5.6 cost crossover", fontsize=4.6, pad=1)
axi.legend(fontsize=3.8, loc="upper left", handlelength=1.2,
           borderpad=0.2, labelspacing=0.25)
axi.text(3.9, 2.2, "L* = alpha/2\n~ 3-4 stages", fontsize=3.8, color="0.3", va="bottom")
axi.set_xlabel("cascade depth L", fontsize=4.2, labelpad=0.5)
axi.set_ylabel("energy (kT D$^2$)", fontsize=4.2, labelpad=0.5)

ax.text(0.995, 0.17,
        "analogue E_A = 2kTL^2D^2\n"
        "hard ceiling 2.35 bits\n"
        "counting E_B = alpha kTLD^2\n"
        "alpha ~ 7, soft ceiling 4.5 bits\n"
        "crossover L* = alpha/2 ~ 3-4",
        fontsize=4.3, ha="right", va="center", color="0.25")
schematic_tag(ax)

# ============ panel b: Theorem T' stratified classification ============
ax = axs[0, 1]
ax.axis("off")
ax.set_xlim(0, 1); ax.set_ylim(0, 1)

ax.text(0.5, 0.985, "Theorem T': stratified classification over degeneracy groups (S5.2)",
        fontsize=6.4, ha="center", va="top", color="0.15", fontweight="bold")

# nested boxes, left half
ax.add_patch(FancyBboxPatch((0.02, 0.06), 0.46, 0.72, boxstyle="round,pad=0.008",
                            fc=OI["verm"], ec=OI["verm"], alpha=0.10, lw=0.9))
ax.add_patch(FancyBboxPatch((0.07, 0.13), 0.36, 0.52, boxstyle="round,pad=0.008",
                            fc=OI["orange"], ec=OI["orange"], alpha=0.12, lw=0.9))
ax.add_patch(FancyBboxPatch((0.12, 0.20), 0.26, 0.32, boxstyle="round,pad=0.008",
                            fc=OI["blue"], ec=OI["blue"], alpha=0.14, lw=0.9))
ax.text(0.045, 0.755, "G2: y(t) -> m(y(t)), monotone reparam",
        fontsize=4.7, color=OI["verm"], fontweight="bold", va="top")
ax.text(0.09, 0.625, "G1: y -> ay + b (add baseline)", fontsize=4.7,
        color=OI["orange"], fontweight="bold", va="top")
ax.text(0.14, 0.495, "G0: y -> ay (multiplicative)", fontsize=4.7,
        color=OI["blue"], fontweight="bold", va="top")
ax.text(0.25, 0.335, "survivors: counts\nproven layer (S5.1)",
        fontsize=4.9, ha="center", color="0.25")
ax.text(0.25, 0.175, "ratio statistics die here;\nlog-linearisation stops",
        fontsize=4.6, ha="center", color="0.4")
ax.text(0.25, 0.085, "survivors: event-time / ordinal\n"
        "max invariant R_y = {(s,t): y(s)<=y(t)} (L1')",
        fontsize=4.4, ha="center", color="0.4")
ax.text(0.25, 0.012, "G0 subset G1 subset G2", fontsize=5.0, ha="center",
        color="0.3", fontweight="bold", va="bottom")

# numeric pins table, right half
rows = [
    ["plateau-free\n(n=300)", "0.053671", "0.053472"],
    ["40% plateau\nnon-tied", "0.0608", "0.0604"],
    ["tied-pair\nflip rate", "1.0000", "1.0000"],
    ["40% plateau\ntotal", "0.1840", "0.1843"],
]
tbl = ax.table(cellText=rows,
               colLabels=["S5.2 pin", "theory", "empirical"],
               bbox=[0.50, 0.54, 0.49, 0.34], cellLoc="center")
tbl.auto_set_font_size(False)
tbl.set_fontsize(4.6)
for (r, c), cell in tbl.get_celld().items():
    cell.set_edgecolor("0.6")
    cell.set_linewidth(0.5)
    if r == 0:
        cell.set_text_props(fontweight="bold")
        cell.set_facecolor("#EEEEEE")

ax.text(0.745, 0.47, "open items (honest list):", fontsize=5.0,
        color=GREY, ha="center", fontweight="bold")
ax.text(0.745, 0.33,
        "noise bound for L1 in general\nfunctional norms remains open\n\n"
        "physical realism of the G2\npremise (limiting narrative)\n\n"
        "event-time jitter bound: stated,\nnot yet written in",
        fontsize=4.8, color=GREY, ha="center", va="center", style="italic")
ax.text(0.70, 0.115, "amplitude unrecoverable from R_y:\n"
        "amplitude must be paid for", fontsize=4.8, ha="center", color="0.3")
schematic_tag(ax)

# ============ panel c: conjecture decision tree ============
ax = axs[1, 0]
ax.axis("off")
ax.set_xlim(0, 1); ax.set_ylim(0, 1)

ax.text(0.5, 0.985, "conjecture decision tree (S1.5 legs, S5.4 conjecture)",
        fontsize=6.4, ha="center", va="top", color="0.15", fontweight="bold")

box(ax, 0.28, 0.80, 0.44, 0.12,
    "system under multiplicative\nscale uncertainty (degeneracy group G)",
    "#DCE9F5", ec=OI["blue"], fs=5.4)

box(ax, 0.22, 0.50, 0.56, 0.22,
    "three-leg criterion for crossing an interface\n"
    "1 algebraic survival: X in T_i\n"
    "2 symbolic reachability: I(X; Y_i) > 0\n"
    "3 temporal bandwidth: tau_X in BW(phi_i)",
    "#EFEFEF", ec="0.4", fs=5.0)
arr(ax, 0.50, 0.80, 0.50, 0.72)

box(ax, 0.02, 0.20, 0.30, 0.19,
    "all legs pass, invariant and\nlikelihood-relevant:\ncounting statistics sufficient",
    "#DDEEDD", ec=OI["green"], fs=5.0)
box(ax, 0.36, 0.20, 0.29, 0.19,
    "leg fails or absolute scale\nreferenced: payment required\n(class III reference, maintained)",
    "#F6E3DE", ec=OI["verm"], fs=5.0)
box(ax, 0.70, 0.20, 0.29, 0.19,
    "third-system test, Morris-Lecar\n(code 25): 2/3 passed, cosine\n0.0647 fails: OPEN branch",
    "#FFFFFF", ec=GREY, tc=GREY, fs=5.0)

arr(ax, 0.36, 0.50, 0.17, 0.39)
arr(ax, 0.50, 0.50, 0.50, 0.39)
arr(ax, 0.64, 0.50, 0.845, 0.39, color=GREY)

ax.text(0.5, 0.125, "conjecture (S5.4): invariant + likelihood-relevant partition\n"
        "statistics must be of counting type; third-system existence test open",
        fontsize=4.9, ha="center", va="top", color="0.3")
schematic_tag(ax)

# ============ panel d: Theorem N noise-robustness sketch ============
ax = axs[1, 1]
t = np.linspace(0.0, 0.42, 400)
for n, c in [(300, OI["blue"]), (1000, OI["verm"])]:
    y = np.minimum(2 * np.exp(-n * t ** 2 / 2), 1.0)
    ax.plot(t, y, "-", lw=1.1, color=c, label=f"n = {n}")
ax.axhline(0.1840, color="0.45", ls="--", lw=0.7)
ax.text(0.415, 0.215, "40% plateau total 0.1840 (empirical 0.1843)",
        fontsize=4.8, color="0.4", ha="right", va="bottom")
ax.set_yscale("log")
ax.set_ylim(1e-6, 3)
ax.set_xlim(0, 0.42)
ax.set_xlabel("deviation t")
ax.set_ylabel("tail bound 2 exp(-n t$^2$/2)  (S5.1)")
ax.legend(fontsize=5.2, loc="upper right")
ax.text(0.98, 0.42,
        "atom wall (N3): tied mass rho gives floor Theta(rho^2),\n"
        "independent of method\n"
        "p53 anchor (s = 15% of amplitude):\n"
        "pulse-baseline flip Phi(-4.71) ~ 1.2e-6",
        transform=ax.transAxes, fontsize=4.9, color="0.3",
        ha="right", va="top")
ax.set_title("Theorem N: concentration is free, ties set the floor", pad=2)
schematic_tag(ax)

for i, ax in enumerate(axs.ravel()):
    panel_label(ax, "abcd"[i])

fig.suptitle("EDFig. 6  Digital-limit theorem diagrams and conjecture decision tree",
             fontsize=8, y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.955])
save(fig, os.path.join(OUT, "EDFig6_digital_limit"))
