import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from edstyle import setup, panel_label, schematic_tag, save, OI
import matplotlib.pyplot as plt

setup()
ROOT = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911"
OUT = os.path.join(ROOT, r"05_主线纲领与设计\论文_细胞通讯审计_2026-09-23\figures_svg\ED")

fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.8))
axs = axs.ravel()

# ---- a: NF discrimination diagram: A^2 linear in (s - s_H) ----
ax = axs[0]
rng = np.random.default_rng(7)
ds = np.linspace(0, 6, 14)
A2 = 0.083 * ds
ax.plot(ds, A2 * (1 + rng.normal(0, 0.03, ds.size)), "o", ms=3, color=OI["blue"], label="NF Hopf (repr.)")
ax.plot(ds, 0.083 * ds, "-", color=OI["blue"], lw=1.0)
ds2 = np.linspace(0, 6, 10)
ax.plot(ds2, 0.055 * ds2 * (1 + rng.normal(0, 0.04, ds2.size)), "s", ms=2.8, color=OI["sky"], label="delayed NF")
ax.plot(ds2, 0.055 * ds2, "--", color=OI["sky"], lw=1.0)
ax.axvline(0, color="0.4", lw=0.7, ls=":")
ax.annotate("s_H = 16.4449", (0.15, 0.02), fontsize=5.5, rotation=90, color="0.35", va="bottom")
ax.text(2.6, 0.10, "R2 = 0.9926 (NF)\nR2 = 0.971 (delayed)", fontsize=5.5, color="0.25")
ax.set_xlabel(r"drive above onset, $s - s_H$"); ax.set_ylabel(r"squared amplitude $A^2$ (au)")
ax.set_title("NF discriminator: A^2 linear in drive")
ax.legend(fontsize=5.2, loc="upper left")
schematic_tag(ax)

# ---- b: amplitude grading: EXC flat vs NF steep vs ML Type-II ----
ax = axs[1]
drive = np.linspace(1, 4.6, 20)
exc = 1.0 * (1 + 0.065 * (drive - 1) / 3.6) * (1 + rng.normal(0, 0.01, drive.size))
nf = 0.28 * np.sqrt(drive) * (1 + rng.normal(0, 0.02, drive.size))
ml = 1.35 - 0.258 * np.log10(drive)
ax.plot(drive, exc, "o-", ms=3, color=OI["green"], label="EXC (FHN): +6.5% over 4.6x, a = 0.014")
ax.plot(drive, nf / nf[0], "s-", ms=3, color=OI["verm"], label="NF: +360% (contrast)")
ax.plot(drive, ml, "^-", ms=3, color=OI["blue"], label="ML Type-II: a = -0.258")
ax.set_xlabel("drive (fold of threshold)"); ax.set_ylabel("amplitude (normalized)")
ax.set_title("amplitude grading separates EXC from NF")
ax.legend(fontsize=5.0, loc="upper left")
schematic_tag(ax)

# ---- c: subcritical normal form: onset + hysteresis, mu sweep ----
ax = axs[2]
def subcrit_branch(r):
    # normal form dx/dt = x(r - 1 + x^2 - x^4)-like amplitude envelope (schematic)
    up = np.sqrt(0.5 + np.sqrt(np.maximum(r - 0.75, 0)))
    lo = np.where(r > 1.0, 0.0, np.nan)
    return lo, up
r = np.linspace(0.74, 1.35, 400)
lo, up = subcrit_branch(r)
ax.plot(r, up, "-", color=OI["blue"], lw=1.2, label="mu = 0 branch")
ax.plot(r[r > 1.0], np.zeros_like(r[r > 1.0]), "-", color=OI["blue"], lw=1.2)
# middle (unstable) branch
mid = np.sqrt(np.maximum(0.5 - np.sqrt(np.maximum(r - 0.75, 0)), 0))
ax.plot(r[r < 1.02], mid[r < 1.02], "--", color=OI["blue"], lw=0.9)
r2 = r + 0.13
lo2, up2 = subcrit_branch(r2)
mid2 = np.sqrt(np.maximum(0.5 - np.sqrt(np.maximum(r2 - 0.75, 0)), 0))
ax.plot(r2, up2, "-", color=OI["orange"], lw=1.2, label="mu = 0.35 branch")
ax.plot(r2[r2 > 1.13], np.zeros_like(r2[r2 > 1.13]), "-", color=OI["orange"], lw=1.2)
ax.plot(r2[r2 < 1.15], mid2[r2 < 1.15], "--", color=OI["orange"], lw=0.9)
ax.annotate("", xy=(1.00, 0.35), xytext=(0.775, 0.35), arrowprops=dict(arrowstyle="<->", lw=0.7, color=OI["verm"]))
ax.text(0.24, 0.42, "hysteresis 0.275 (predicted 0.25)", fontsize=5, color=OI["verm"], transform=ax.transAxes)
ax.annotate("", xy=(1.13, 0.55), xytext=(1.00, 0.55), arrowprops=dict(arrowstyle="->", lw=0.7, color="0.35"))
ax.text(0.66, 0.26, "onset 1.00 -> 1.13\n(mu: 0 -> 0.35)", fontsize=5, color="0.35", transform=ax.transAxes)
ax.set_xlabel("normal-form parameter r"); ax.set_ylabel("oscillation amplitude (au)")
ax.set_title("subcritical normal form: onset and hysteresis")
ax.legend(fontsize=5.2, loc="upper left")
schematic_tag(ax)

# ---- d: counting law ----
ax = axs[3]
D = np.linspace(0.05, 0.8, 24)
Dc = 0.261
T = 5.5 * np.sqrt(np.maximum(D - Dc, 1e-3) / (1 - Dc)) * 2.14 / 2.14
Traw = 2.14 * np.sqrt(np.maximum(D - Dc, 1e-3)) * 5.5
y = Traw / Traw.max() * 2.30
ax.plot(D, y * (1 + rng.normal(0, 0.03, D.size)), "o", ms=3, color=OI["blue"])
ax.plot(D, y, "-", color=OI["blue"], lw=1.0)
ax.axvline(Dc, color="0.4", ls=":", lw=0.7)
ax.text(0.285, 1.9, "D_c = 0.261", fontsize=5.5, rotation=90, color="0.35")
ax.text(0.42, 0.35, "R2 = 0.920\nT inversion 2.14 au\nmeasured 2.30 au", fontsize=5.5, color="0.25")
ax.set_xlabel("damage dose D (au)"); ax.set_ylabel("mean inter-pulse interval T (au)")
ax.set_title("counting law above damage threshold")
schematic_tag(ax)

# ---- e: Batchelor two arms ----
ax = axs[4]
ncs = np.array([137.9, 193.1, 137.9])
uv = np.array([84.8, 103.3, 181.5, 366.3, 331.5])
xg1 = np.arange(ncs.size); xg2 = np.arange(uv.size) + ncs.size + 1
ax.bar(xg1, ncs, width=0.7, color=OI["blue"], label="NCS (1.40x)")
ax.bar(xg2, uv, width=0.7, color=OI["verm"], label="UV (4.3x)")
ax.set_xticks(list(xg1) + list(xg2))
ax.set_xticklabels(["N1", "N2", "N3", "U1", "U2", "U3", "U4", "U5"], fontsize=5.5)
ax.set_ylabel("p53 pulse amplitude (AU)")
ax.set_title("Batchelor two arms: amplitude invariance vs grading")
ax.text(0.02, 0.96, "a_NCS = -0.001 +/- 0.113\na_UV  = +0.981 +/- 0.114\n"
        "separation 6.1 sigma (independent errors)\n8.6 sigma (frozen-card one-sided)",
        transform=ax.transAxes, fontsize=5.2, va="top", color="0.25")
ax.set_ylim(0, 430)
ax.legend(fontsize=5.2, loc="upper right")
schematic_tag(ax)

# ---- f: TrackC2 signature grid (alpha, nu) ----
ax = axs[5]
pts = [("gamma-IR (DSB)", 0.05, 0.10, 1.4, 0.4, OI["green"]),
       ("NCS (DSB)", -0.0008, 0.113, 2.0, 0.5, OI["blue"]),
       ("UV", 0.981, 0.114, 0.0, 0.2, OI["verm"])]
for name, a, ae, nu, ne, c in pts:
    ax.errorbar(a, nu, xerr=ae, yerr=ne, fmt="o", ms=4.5, color=c, capsize=2.5, lw=0.9, label=name)
ax.plot(1.135, 0.0, "s", ms=5, color="0.2", mfc="none", label="model anchor NF")
ax.plot(0.04, 1.4, "D", ms=4.5, color=OI["pink"], mfc="none", label="model anchor EXC (nu = tau_r/T)")
ax.axhline(0, color="0.5", lw=0.6); ax.axvline(0, color="0.5", lw=0.6)
ax.text(0.55, 2.3, "counting regime\n(a ~ 0, nu > 0)", fontsize=5.5, color="0.4", ha="center")
ax.text(0.30, 0.07, "amplitude regime\n(a ~ 1, nu ~ 0)", fontsize=5.5, color="0.4", ha="center", transform=ax.transAxes)
ax.text(1.42, 0.50, "tau_r ~ 7.7 h (= nu*T)", fontsize=5, color="0.4", ha="right")
ax.set_xlabel(r"amplitude-dose exponent $\alpha$"); ax.set_ylabel(r"counting exponent $\nu$")
ax.set_xlim(-0.45, 1.45); ax.set_ylim(-0.75, 2.75)
ax.set_title("TrackC2 damage signature grid (code 52)")
ax.legend(fontsize=4.8, loc="center right")
schematic_tag(ax)

for i, ax in enumerate(axs):
    panel_label(ax, "abcdef"[i])

fig.suptitle("EDFig. 3  p53 discrimination diagram, normal-form hysteresis, counting law, "
             "Batchelor two-arm audit,\nTrackC2 damage signature grid (schematic redraws of archived results)",
             fontsize=8)
fig.tight_layout(rect=[0, 0, 1, 0.925])
save(fig, os.path.join(OUT, "EDFig3_p53_discrimination_grid"))
