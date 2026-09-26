import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from edstyle import setup, panel_label, schematic_tag, save, OI
import matplotlib.pyplot as plt

setup()
ROOT = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911"
OUT = os.path.join(ROOT, r"05_主线纲领与设计\论文_细胞通讯审计_2026-09-23\figures_svg\ED")
rng = np.random.default_rng(11)

fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.8))
axs = axs.ravel()

# ---- a: soft regime: A^2 linear in C, Hopf at C_H ----
ax = axs[0]
C = np.linspace(2.49e-3, 8e-3, 14)
A2 = 0.35 * (C - 2.49e-3) / 1e-3
ax.plot(C * 1e3, A2 * (1 + rng.normal(0, 0.02, C.size)), "o", ms=3, color=OI["blue"])
ax.plot(C * 1e3, A2, "-", color=OI["blue"], lw=1.0)
ax.axvline(2.49, color="0.4", ls=":", lw=0.7)
ax.text(2.6, 0.08, "C_H = 2.49e-3", fontsize=5.5, rotation=90, color="0.35", va="bottom")
ax.text(4.4, 0.30, "R2 = 0.996\na = +0.729\nno hysteresis", fontsize=5.5, color="0.25")
ax.set_xlabel(r"IKK drive $C$ ($\times 10^{-3}$)"); ax.set_ylabel(r"squared amplitude $A^2$ (au)")
ax.set_title("soft regime: supercritical Hopf (B = 3, d = 0.005)")
schematic_tag(ax)

# ---- b: spiky regime bifurcation: double Hopf + canard ----
ax = axs[1]
Cx = np.logspace(np.log10(4e-3), np.log10(0.25), 300)
amp = np.where(Cx < 5.89e-3, 0.02,
       np.where(Cx < 6.00e-3, 0.061,
       np.where(Cx < 0.16741, 0.175 * (Cx / 6e-3) ** 0.35, 0.02)))
amp = np.clip(amp, 0.02, None)
ax.plot(Cx, amp, "-", color=OI["verm"], lw=1.2)
for xh, lab in [(5.89e-3, "H1 0.00589"), (0.16741, "H2 0.16741")]:
    ax.axvline(xh, color="0.4", ls=":", lw=0.7)
    ax.text(xh * 1.06, 0.018, lab, fontsize=5, rotation=90, color="0.35")
ax.annotate("", xy=(6.0e-3, 0.175), xytext=(6.0e-3, 0.061),
            arrowprops=dict(arrowstyle="->", lw=0.9, color=OI["blue"]))
ax.text(0.24, 0.11, "canard at C ~ 0.00600\nA: 0.061 -> 0.175", fontsize=5, color=OI["blue"], transform=ax.transAxes)
ax.text(0.02, 0.97, "110-pt bidirectional continuation:\nmax fwd-rev diff 2.05e-3 (zero hysteresis)",
        transform=ax.transAxes, fontsize=5, va="top", color="0.25")
ax.set_xscale("log")
ax.set_xlabel("IKK drive C"); ax.set_ylabel("oscillation amplitude (au)")
ax.set_title("spiky regime: double Hopf + canard (B = 954.5)")
schematic_tag(ax)

# ---- c: spiky regime exponent decomposition + period/duty ----
ax = axs[2]
seg1 = np.linspace(6e-3, 0.02, 8)     # rising
seg2 = np.linspace(0.02, 0.06, 6)     # plateau-ish peak
seg3 = np.linspace(0.06, 0.167, 8)    # falling
a1 = 0.175 * (seg1 / 6e-3) ** 0.549
a3 = 0.175 * (seg3 / 0.167) ** 0.438 * 1.0
ax.plot(seg1, a1, "o-", ms=3, color=OI["green"])
ax.plot(seg3, a3 * (0.175 / a3[-1]), "s-", ms=3, color=OI["verm"])
ax.plot([6e-3, 0.10], [0.061, 0.24], "--", color="0.5", lw=0.8)
ax.text(0.0066, 0.255, "rising: a = +0.549", fontsize=5.2, color=OI["green"], rotation=24)
ax.text(0.052, 0.075, "falling: a = -0.438", fontsize=5.2, color=OI["verm"])
ax.text(0.016, 0.078, "incl. canard: a = +0.844", fontsize=5.2, color="0.45", rotation=20)
ax.set_xscale("log")
ax.set_xlabel("IKK drive C"); ax.set_ylabel("amplitude (au)")
ax.text(0.97, 0.93, "period 0.37-3.20 (8.8x)\nduty cycle 10.8% at C = 0.035",
        transform=ax.transAxes, fontsize=5.2, va="top", ha="right", color="0.25")
ax.set_title("spiky regime: exponent sign flips at peak")
schematic_tag(ax)

# ---- d: Tay 2010 benchmark ----
ax = axs[3]
cats = ["amplitude\nexponent a", "counting slope nu\n(per decade)"]
model = [0.729, 0.0]
tay = [0.151, 0.58]
x = np.arange(2); w = 0.32
ax.bar(x - w / 2, model, w, color=OI["blue"], label="model (spiky regime)")
ax.bar(x + w / 2, tay, w, color=OI["orange"], label="Tay 2010 experiment")
ax.axhspan(0, 0.83, xmin=0.5, xmax=1.0, color=OI["green"], alpha=0.10)
ax.text(1.0, 0.90, "nu upper bound 0.83", fontsize=5, ha="center", color=OI["green"])
for xi, v in zip(x - w / 2, model):
    ax.text(xi, v + 0.10 if v == 0 else v + 0.02, f"{v:.3f}", ha="center", fontsize=5.5)
for xi, v in zip(x + w / 2, tay):
    ax.text(xi, v + 0.02, f"{v:.3f}", ha="center", fontsize=5.5)
ax.set_xticks(x); ax.set_xticklabels(cats, fontsize=6)
ax.set_ylabel("value")
ax.set_title("model vs Tay 2010 benchmark")
ax.legend(fontsize=5.2, loc="upper left")
ax.set_ylim(0, 1.05)
schematic_tag(ax)

# ---- e: code 12 patch audit (rejected) ----
ax = axs[4]
names = ["fold change\n(4 sets)", "soft regime\noscillation", "spiky window\na", "amplification\n(resp/pert)"]
vals = [0, 0, -0.005, 1.0]
cols = [OI["verm"], OI["verm"], OI["verm"], OI["verm"]]
b = ax.bar(np.arange(4), [1, 1, 1, 1], color="0.92", width=0.62)
labels = ["0 / 150", "abolished", "-0.005", "3000 / 3000 = 1"]
for i, t in enumerate(labels):
    ax.text(i, 0.5, t, ha="center", va="center", fontsize=5.8, color=OI["verm"])
ax.set_xticks(np.arange(4)); ax.set_xticklabels(names, fontsize=5.2)
ax.set_yticks([])
ax.set_title("code 12 patch: all four checks fail")
ax.text(0.02, 0.13, "spiky window moved to (9.11e-3, 0.20) but exponent lost", transform=ax.transAxes,
        fontsize=5, color="0.4")
schematic_tag(ax)

# ---- f: code 15 patch audit (rejected): nu vs K_T ----
ax = axs[5]
KT = np.logspace(np.log10(0.005), np.log10(0.2), 40)
for cmax, col in zip([0.004, 0.008, 0.02], [OI["sky"], OI["blue"], OI["pink"]]):
    nu = 1.6 - 2.6 * cmax - 0.55 * np.log10(KT / 0.005) * (0.4 + 12 * cmax)
    ax.plot(KT, nu, "-", color=col, lw=1.1, label=f"C_max = {cmax:g}")
ax.axhspan(0, 0.83, color=OI["green"], alpha=0.10)
ax.text(0.006, 0.62, "Tay band: nu <= 0.83", fontsize=5, color=OI["green"])
ax.plot(0.03, 0.900, "o", ms=5, color=OI["verm"])
ax.annotate("min nu = 0.900\n(C_max = 0.008, K_T = 0.03)", (0.03, 0.900),
            xytext=(0.0053, 1.00), fontsize=5, color=OI["verm"],
            arrowprops=dict(arrowstyle="->", lw=0.7, color=OI["verm"]))
ax.set_xscale("log")
ax.set_xlabel(r"$K_T$"); ax.set_ylabel(r"counting slope $\nu$ (per decade)")
ax.set_title("code 15: 8/24 pass a, none pass nu")
ax.legend(fontsize=5.2, loc="center right", bbox_to_anchor=(1.0, 0.33))
schematic_tag(ax)

for i, ax in enumerate(axs):
    panel_label(ax, "abcdef"[i])

fig.suptitle("EDFig. 7  NF-kB two-regime audit and the two rejected patches "
             "(codes 12, 15; schematic redraws of archived results)", fontsize=8)
fig.tight_layout(rect=[0, 0, 1, 0.94])
save(fig, os.path.join(OUT, "EDFig7_NFkB_two_regime_audit"))
