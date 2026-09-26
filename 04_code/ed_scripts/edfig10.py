import json, sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from edstyle import setup, panel_label, save, OI
import matplotlib.pyplot as plt

setup()
ROOT = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911"
OUT = os.path.join(ROOT, r"05_主线纲领与设计\论文_细胞通讯审计_2026-09-23\figures_svg\ED")
D = json.load(open(os.path.join(ROOT, r"05_主线纲领与设计\论文_细胞通讯审计_2026-09-23\figures_svg\ED\scripts\chemo54_rerun.json")))

fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.6))
axs = axs.ravel()

# ---- a: Pareto sweep + same-source control ----
ax = axs[0]
fr = D["front"]
ax.add_patch(plt.Rectangle((0.9, 0), 0.12, 0.1, color=OI["green"], alpha=0.15, zorder=0))
ax.text(0.905, 0.02, "acceptance\nregion\nR2>0.9, err<0.1 dex", fontsize=5, color=OI["green"], va="bottom")
ax.plot([f["r2"] for f in fr], [f["dk"] for f in fr], "o-", color=OI["blue"], ms=4.5, zorder=3)
for f in fr:
    ax.annotate(f"lam={f['lam']:g}", (f["r2"], f["dk"]), textcoords="offset points",
                xytext=(5, 4), fontsize=5, color=OI["blue"])
# same-source (code 62) and whole-population control, SI S8.5
ss = [(0.971, 0.573), (0.868, 0.205), (0.599, 0.071)]
wp = [(0.976, 0.624), (0.848, 0.261), (0.525, 0.074)]
ax.plot(*zip(*ss), "s--", color=OI["verm"], ms=3.5, mfc="none", label="same-source (code 62)")
ax.plot(*zip(*wp), "d:", color=OI["grey"], ms=3.5, mfc="none", label="whole-population control")
ax.set_xlabel(r"amplitude $R^2$"); ax.set_ylabel(r"K$_{1/2}$ mean |$\Delta$log$_{10}$| (dex)")
ax.set_title("Pareto front: no point enters the acceptance region")
ax.legend(fontsize=5.2, loc="upper left")
ax.set_xlim(0.45, 1.05); ax.set_ylim(0, 0.7)

# ---- b: repair-path grid ----
ax = axs[1]
ax.add_patch(plt.Rectangle((0.9, 0), 0.12, 0.1, color=OI["green"], alpha=0.15, zorder=0))
quant = [(2, 0.953, 0.182, "M2c imperfect adaptation"),
         (4, 0.950, 0.257, "V3 readout nonlinearity"),
         (6, 0.680, 0.052, "gain retuning per background"),
         (7, 0.946, 0.179, "code 60 heterogeneous population")]
for n, r2, dk, lab in quant:
    ax.plot(r2, dk, "o", ms=6, color=OI["pink"])
    ax.annotate(f"#{n}", (r2, dk), textcoords="offset points", xytext=(6, -2), fontsize=6,
                color=OI["pink"], fontweight="bold")
ax.plot(fr[0]["r2"], fr[0]["dk"], "o", ms=6, color=OI["blue"])
ax.annotate("standard (lam=0)", (fr[0]["r2"], fr[0]["dk"]), textcoords="offset points",
            xytext=(-4, 7), fontsize=5, color=OI["blue"], ha="right")
ax.plot(fr[-1]["r2"], fr[-1]["dk"], "o", ms=6, color=OI["blue"])
ax.annotate("standard (lam=30)", (fr[-1]["r2"], fr[-1]["dk"]), textcoords="offset points",
            xytext=(6, -2), fontsize=5, color=OI["blue"])
notes = ("#1 two subpopulations: directional failure (alpha -> 0.05-0.08)\n"
         "#3 protocol mismatch: excluded as primary cause\n"
         "#5 TCS depletion: front digit-for-digit unchanged\n"
         "#8 asymmetric N: n(B) 1->17 within single-N envelope")
ax.text(0.03, 0.60, notes, transform=ax.transAxes, fontsize=5, va="top", color="0.25")
ax.set_xlabel(r"amplitude $R^2$"); ax.set_ylabel(r"K$_{1/2}$ error (dex)")
ax.set_title("eight repair paths: all outside the acceptance region")
ax.set_xlim(0.45, 1.05); ax.set_ylim(0, 0.7)

# ---- c: Weber line + residuals ----
ax = axs[2]
B = np.array(D["K12_B"]); med = np.array(D["K12_obs"]); iqr = np.array(D["K12_iqr"])
ax.errorbar(B, med, yerr=np.vstack([med - iqr[:, 0], iqr[:, 1] - med]), fmt="o",
            color=OI["blue"], ms=4, elinewidth=0.6, capsize=1.5, label="measured K$_{1/2}$ (IQR)")
Bb = np.logspace(-2.3, 2.1, 100)
ax.plot(Bb, 1.17 * (1.95 + Bb), "--", color=OI["verm"], lw=0.9, label=r"Weber line 1.17(1.95 + $B$)")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"background $B$ (uM MeAsp)"); ax.set_ylabel(r"K$_{1/2}$ (uM)")
ax.set_title("Weber line: max deviation 0.10 dex, mean 0.03 dex")
ax.legend(fontsize=5.2, loc="upper left")
axin = ax.inset_axes([0.42, 0.08, 0.55, 0.38])
res = np.log10(med / (1.17 * (1.95 + B)))
axin.bar(range(7), res, color=OI["sky"], width=0.7)
axin.axhline(0, color="0.4", lw=0.5)
axin.set_xticks(range(7), [""]*7)
axin.tick_params(labelsize=4)
axin.set_xlabel("B (uM)", fontsize=4.5)
axin.set_ylabel("residual (dex)", fontsize=5)
for sp in axin.spines.values():
    sp.set_linewidth(0.4)

# ---- d: censoring per background ----
ax = axs[3]
cens = np.array([D["cens"][str(float(b))] for b in B])
tot = np.array([D["tot"][str(float(b))] for b in B])
est = tot - cens
xx = np.arange(7)
ax.bar(xx, tot, color="0.8", width=0.65, label="right-censored")
ax.bar(xx, est, color=OI["blue"], width=0.65, label="estimable")
for i, (e, t_) in enumerate(zip(est, tot)):
    ax.text(i, t_ + 4, f"{100*cens[i]/t_:.0f}%", ha="center", fontsize=5, color="0.4")
ax.set_xticks(xx, [f"{b:g}" for b in B])
ax.set_xlabel(r"background $B$ (uM)"); ax.set_ylabel("cells")
ax.set_title("right-censoring: 45-68% of cells per background")
ax.legend(fontsize=5.2, loc="upper left")

# ---- e: prospective midpoint clauses ----
ax = axs[4]
pred = {"V1": 388.0, "V2": 122.0}
meas = {"V1": (401.3, 388.5, 407.6), "V2": (133.2, 122.8, 140.4), "V2_alt": (130.8, None, None)}
ax.plot([30, 600], [30, 600], color="0.4", lw=0.6)
for tol in [0.15]:
    ax.fill_between([30, 600], [30 * 10**-tol, 30 * 10**-tol], [600 * 10**tol, 600 * 10**tol],
                    color=OI["green"], alpha=0.12)
ax.text(560, 70, "+/-0.15 dex frozen tolerance", fontsize=5, color=OI["green"], ha="right", va="bottom")
cols = {"V1": OI["blue"], "V2": OI["verm"], "V2_alt": OI["grey"]}
for k, (m, lo, hi) in meas.items():
    p = pred.get(k, 122.0)
    yerr = None if lo is None else [[m - lo], [hi - m]]
    ax.errorbar(p, m, yerr=yerr, fmt="o", color=cols[k], ms=5, capsize=2, elinewidth=0.7,
                label=f"{k}: pred {p:g}, meas {m:g}")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xticks([100, 200, 400]); ax.set_yticks([100, 200, 400])
from matplotlib.ticker import ScalarFormatter
ax.xaxis.set_major_formatter(ScalarFormatter()); ax.yaxis.set_major_formatter(ScalarFormatter())
ax.minorticks_off()
ax.set_xlabel(r"predicted K$_{1/2}$ (uM, effective axis)"); ax.set_ylabel(r"measured K$_{1/2}$ (uM)")
ax.set_title("prospective arms (code 70): P-2 hit 2/2")
ax.legend(fontsize=5.2, loc="upper left")
ax.set_xlim(60, 600); ax.set_ylim(60, 600)

# ---- f: prospective amplitude + structural clauses ----
ax = axs[5]
xx = np.arange(3)
amp_dev = [0.105, 0.126, 0.069]
str_dev = [-0.005, 0.177, 0.202]
b1 = ax.bar(xx - 0.18, amp_dev, width=0.36, color=OI["orange"], label="P-1b max amplitude dev")
b2 = ax.bar(xx + 0.18, str_dev, width=0.36, color=OI["sky"], label=r"P-1a $R^2$(logF) - $R^2$(logr)")
ax.axhline(0.08, color=OI["orange"], ls="--", lw=0.7)
ax.text(-0.42, 0.085, "amp tolerance 0.08", fontsize=5, color=OI["orange"])
ax.axhline(0.30, color=OI["blue"], ls="--", lw=0.7)
ax.text(-0.42, 0.305, "structural line +0.3", fontsize=5, color=OI["blue"])
for i, v in enumerate(amp_dev):
    ax.text(i - 0.18, v + 0.006, f"{v:.3f}", ha="center", fontsize=5)
for i, v in enumerate(str_dev):
    ax.text(i + 0.18, v + 0.006, f"{v:+.3f}", ha="center", fontsize=5)
ax.set_xticks(xx, ["V1", "V2", "V2_alt (control)"])
ax.set_ylabel("deviation")
ax.set_title("P-1b missed at low F; P-1a underpowered")
ax.legend(fontsize=5.2, loc="upper right")
ax.set_ylim(-0.03, 0.36)

for i, ax in enumerate(axs):
    panel_label(ax, "abcdef"[i])

fig.suptitle("EDFig. 10  Chemotaxis receptor-layer rupture: Pareto sweep, repair grid, Weber line, "
             "prospective arms (Moore et al. 2024; codes 54, 62, 70)", fontsize=8)
fig.tight_layout(rect=[0, 0, 1, 0.955])
save(fig, os.path.join(OUT, "EDFig10_chemotaxis_rupture"))
