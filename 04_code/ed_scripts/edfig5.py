import json, sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from edstyle import setup, panel_label, save, OI
import matplotlib.pyplot as plt

setup()
ROOT = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911"
OUT = os.path.join(ROOT, r"05_主线纲领与设计\论文_细胞通讯审计_2026-09-23\figures_svg\ED")
dA = json.load(open(os.path.join(ROOT, r"03_细胞线3\结果\代码51r4A_噪声标定与稳健性.json"), encoding="utf-8"))
dB = json.load(open(os.path.join(ROOT, r"03_细胞线3\结果\代码51r4B_Fisher_M4000.json"), encoding="utf-8"))
d52 = json.load(open(os.path.join(ROOT, r"03_细胞线3\结果\代码52_TrackC2_结果.json"), encoding="utf-8"))

fig = plt.figure(figsize=(7.2, 6.4))
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.15])
axa = fig.add_subplot(gs[0, 0]); axb = fig.add_subplot(gs[0, 1])
axc = fig.add_subplot(gs[1, 0]); axd = fig.add_subplot(gs[1, 1])
axe = fig.add_subplot(gs[2, :])

# ---- a: fast-noise calibration ----
cal = dA["calib"]
sig = [r["sigma"] for r in cal]
axa.plot(sig, [r["ipi_cv"] for r in cal], "o-", color=OI["blue"], ms=3.5, label="intra-cellular IPI CV")
axa.plot(sig, [r["cv_T"] for r in cal], "s--", color=OI["orange"], ms=3.5, label="population CV(T)")
axa.axhline(0.30, color=OI["verm"], ls=":", lw=0.8, label="literature IPI CV = 0.30")
axa.axvline(0.25, color="0.4", ls="--", lw=0.7)
axa.text(0.26, 0.62, r"$\sigma^*$ = 0.25" + "\n(IPI CV = 0.316, N = 10.4)", fontsize=5.5)
axa.set_xscale("log")
axa.set_xlabel(r"fast noise amplitude $\sigma$"); axa.set_ylabel("CV")
axa.set_title("fast-noise calibration")
axa.legend(fontsize=5.2, loc="upper left")
axa.set_ylim(0, 0.75)
ax2 = axa.twinx()
ax2.plot(sig, [r["N_mean"] for r in cal], "^:", color=OI["grey"], ms=3, label="N (right)")
ax2.set_ylabel("pulse count N", color=OI["grey"])
ax2.tick_params(axis="y", labelcolor=OI["grey"], labelsize=6)
ax2.spines["right"].set_linewidth(0.6)

# ---- b: P2' quantitative reproduction at realistic noise ----
lm = dA["lit_matched"]
xm = [r["mu_d_over_T"] for r in lm]; ym = [r["ratio_lit"] for r in lm]
axb.plot(xm, ym, "o-", color=OI["pink"], ms=4, label=r"model at $\sigma^*$ = 0.25")
axb.axhline(5.8, color=OI["verm"], ls=":", lw=0.8, label="literature 5.8 (Lahav)")
axb.axhline(4, color="0.4", ls="--", lw=0.7, label="P2' criterion = 4")
axb.annotate("", xy=(1.0, 5.7), xytext=(0.73, 3.01),
             arrowprops=dict(arrowstyle="->", color="0.3", lw=0.8))
axb.text(0.62, 4.6, r"$\mu_d/T$ = 1.0 extrapolation $\approx$ 5.7", fontsize=5.5)
axb.set_xlabel(r"$\mu_d$ / T"); axb.set_ylabel(r"Var($t_1$) / mean Var(IPI)")
axb.set_title("sensing-stage signature reproduces 5.7 vs measured 5.8")
axb.legend(fontsize=5.2, loc="upper left")
axb.set_ylim(-0.3, 6.3)

# ---- c: M=4000 Fisher spectrum ----
ev = np.array(dB["evals"])
axc.bar(np.arange(1, 8), np.log10(np.maximum(ev, 1e-16)), color=OI["blue"])
axc.set_xlabel("mode index"); axc.set_ylabel(r"log$_{10}$ eigenvalue")
axc.set_title(r"M = 4000 central spectrum: $\lambda_7/\lambda_1 \leq 4.7\times10^{-19}$")

# ---- d: noise-colour plane ----
het = dA["hetero"]
mk = {"a": "o", "eps": "s", "s": "^", "tau_r": "D"}
cc = {"a": OI["blue"], "eps": OI["orange"], "s": OI["green"], "tau_r": OI["pink"]}
for p in ["a", "eps", "s", "tau_r"]:
    rows = [r for r in het if r["param"] == p]
    axd.plot([r["cv_T"] for r in rows], [r["ipi_cv"] for r in rows], mk[p] + "-",
             color=cc[p], ms=3.5, lw=0.8, label=f"slow het. {p}")
cal0 = [r for r in cal if r["sigma"] <= 0.5]
axd.plot([r["cv_T"] for r in cal0], [r["ipi_cv"] for r in cal0], "o--",
         color=OI["verm"], ms=3.5, lw=0.8, label=r"fast white noise $\sigma$")
axd.axhline(0.30, color=OI["verm"], ls=":", lw=0.8)
axd.text(0.155, 0.315, "literature IPI CV = 0.30", fontsize=5.5, color=OI["verm"])
axd.set_xlabel("population CV(T)"); axd.set_ylabel("intra-cellular IPI CV")
axd.set_title("noise-colour plane: slow heterogeneity cannot reach 0.30")
axd.legend(fontsize=5.0, loc="upper left", ncol=2)
axd.set_xlim(0, 0.21); axd.set_ylim(0, 0.48)

# ---- e: TrackC2 discrimination plane ----
axe.axvspan(-0.35, 0.1, color=OI["green"], alpha=0.10)
axe.axvspan(0.5, 1.5, color=OI["verm"], alpha=0.08)
axe.text(-0.33, 7.6, "EXC region (digital)\n" + r"$\alpha \leq 0.1$, $\nu > 0$", fontsize=6, color=OI["green"])
axe.text(1.05, 9.7, "NF region (analogue)\n" + r"$\alpha \geq 0.5$", fontsize=6, color=OI["verm"])
pts = d52["points"]
names = {"γ-IR（DSB）": "gamma-IR (DSB)", "NCS（DSB）": "NCS (DSB)", "UV": "UV"}
cols = {"gamma-IR (DSB)": OI["blue"], "NCS (DSB)": OI["pink"], "UV": OI["orange"]}
for p in pts:
    nm = names[p["name"]]
    lo = 0.0 if nm == "UV" else p["nerr"]
    axe.errorbar(p["alpha"], p["nu"], xerr=p["aerr"],
                 yerr=np.array([[lo], [p["nerr"]]]),
                 fmt="o", color=cols[nm], ms=5, capsize=2, elinewidth=0.7,
                 label=f"{nm}: alpha = {p['alpha']:+.2f} +/- {p['aerr']:.2f}, nu = {p['nu']:.1f}")
for m, mkk, cc2, lab in [(d52["models"][0], "s", OI["verm"], "model anchor NF: alpha = 1.13, nu ~ 0"),
                          (d52["models"][1], "^", OI["green"], "model anchor EXC: alpha = 0.04, nu = tau_r/T")]:
    nu = m["nu"] if m["nu"] is not None else 10.0
    axe.plot(m["alpha"], nu, mkk, color=cc2, ms=6, mfc="none", mew=1.2, label=lab)
axe.set_xlabel(r"$\alpha = \partial \ln A\,/\,\partial \ln D$  (amplitude index)")
axe.set_ylabel(r"$\nu = \partial N\,/\,\partial \ln D$  (counting index)")
axe.set_title("TrackC2: measured (alpha, nu) signatures on the discrimination diagram (code 52)")
axe.legend(fontsize=5.2, loc="center left", ncol=2)
axe.set_xlim(-0.35, 1.5); axe.set_ylim(-1.2, 11)

for ax, L in zip([axa, axb, axc, axd, axe], "abcde"):
    panel_label(ax, L)

fig.suptitle("EDFig. 5  Dispersion signatures and the noise-colour plane (code 51r4; code 52 TrackC2)", fontsize=8)
fig.tight_layout(rect=[0, 0, 1, 0.965])
save(fig, os.path.join(OUT, "EDFig5_dispersion_noise_colour_51r4_52"))
