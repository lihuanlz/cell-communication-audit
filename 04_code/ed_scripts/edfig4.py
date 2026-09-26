import json, sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from edstyle import setup, panel_label, save, OI
import matplotlib.pyplot as plt

setup()
ROOT = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911"
OUT = os.path.join(ROOT, r"05_主线纲领与设计\论文_细胞通讯审计_2026-09-23\figures_svg\ED")
d51 = json.load(open(os.path.join(ROOT, r"03_细胞线3\结果\代码51r3_p53_Fisher审计_结果.json"), encoding="utf-8"))
d4b = json.load(open(os.path.join(ROOT, r"03_细胞线3\结果\代码51r4B_Fisher_M4000.json"), encoding="utf-8"))

# ---------------- panel a: FHN pulse train re-simulation (locked parameter set, SI S4.1) ----------------
th = d51["theta0"]
s, tau_r, a_, b_, eps, I0, sig = th["s"], th["tau_r"], th["a"], th["b"], th["eps"], th["I0"], th["sigma"]
dt = 0.002
T_tot = 50.0
n = int(T_tot / dt)
t = np.arange(n) * dt
rng = np.random.default_rng(51)
v = np.zeros(n); w = np.zeros(n)
v[0], w[0] = -1.0, -0.6
for i in range(n - 1):
    I = I0 * np.exp(-t[i] / tau_r)
    vn, wn = v[i], w[i]
    v[i+1] = vn + dt * s * (vn - vn**3/3 - wn) + sig * np.sqrt(dt) * rng.standard_normal()
    w[i+1] = wn + dt * s * eps * (vn + a_ - b_ * wn - I)
# count pulses (upward crossings of v=0)
cross = np.where((v[:-1] < 0) & (v[1:] >= 0))[0]
Npulse = len(cross)

fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.4))
axs = axs.ravel()

ax = axs[0]
ax.plot(t, v, color=OI["blue"], lw=0.5)
ax.axhline(0, color="0.6", lw=0.4, ls=":")
ax.set_xlabel("time (au)"); ax.set_ylabel("v (FHN fast variable)")
ax.set_title(f"FHN pulse train, synchronous integration, N = {Npulse}")

# ---------------- panel b: Fisher spectra ----------------
ax = axs[1]
ev2 = np.array(d51["fisher_h2"]["evals"])
ev5 = np.array(d51["fisher_h5"]["evals"])
ev4 = np.array(d4b["evals"])
x = np.arange(1, 8)
w_ = 0.27
ax.bar(x - w_, np.log10(np.maximum(ev2, 1e-16)), width=w_, color=OI["blue"], label="M=400 central, h=2%")
ax.bar(x, np.log10(np.maximum(ev5, 1e-16)), width=w_, color=OI["sky"], label="M=400 central, h=5%")
ax.bar(x + w_, np.log10(np.maximum(ev4, 1e-16)), width=w_, color=OI["verm"], label="M=4000 central, h=2%")
ax.set_xlabel("mode index"); ax.set_ylabel(r"log$_{10}$ eigenvalue")
ax.set_title(r"Fisher spectrum: $\lambda_7/\lambda_1 \leq 4.7\times10^{-19}$")
ax.legend(loc="lower left", fontsize=5.5)

# ---------------- panel c: eigenvector heatmap (h2) ----------------
ax = axs[2]
ev = np.abs(np.array(d51["fisher_h2"]["evecs"]))
im = ax.imshow(ev, cmap="viridis", aspect="auto", vmin=0, vmax=1)
ax.set_xticks(range(7), [f"v{i}" for i in range(1, 8)])
ax.set_yticks(range(7), d51["param_names"])
ax.set_title("eigenvectors (softest v7: sigma 0.996)")
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
cb.ax.tick_params(labelsize=5.5)
cb.outline.set_linewidth(0.6)

# ---------------- panel d: channel shares ----------------
ax = axs[3]
share = np.array(d51["fisher_h2"]["share"]) * 100
obs = d51["obs_names"]
cols = [OI["green"]] * 6
cols[2] = OI["verm"]  # amplitude highlighted
bars = ax.bar(obs, share, color=cols)
for b, sv in zip(bars, share):
    ax.text(b.get_x() + b.get_width()/2, sv + 1, f"{sv:.1f}", ha="center", fontsize=6)
ax.set_ylabel("Fisher information share (%)")
ax.set_title("channel shares: timing-type 94.2%, amplitude 5.5%")
ax.set_ylim(0, 70)

# ---------------- panel e: upstream delay scan ----------------
ax = axs[4]
rows = d51["delay"]["rows"]
mud = [r["mu_d_over_T"] for r in rows]
vr = [r["var_ratio"] for r in rows]
ax.plot(mud, vr, "o-", color=OI["pink"], ms=4)
ax.axhline(4, color="0.4", ls="--", lw=0.7, label="P2' criterion = 4")
ax.axvline(0.73, color=OI["verm"], ls=":", lw=0.7)
ax.text(0.77, 20, "literature calibration\n$\\mu_d/T$ = 0.73 (240/330 min),\nratio = 1953", fontsize=5.5, color=OI["verm"], va="center")
ax.set_yscale("log")
ax.set_xlabel(r"$\mu_d$ / T"); ax.set_ylabel(r"Var($t_1$) / Var(T)")
ax.set_title("stochastic sensing stage flips the dispersion ratio")
ax.legend(loc="lower right", fontsize=5.5)

# ---------------- panel f: counting law + gain dispersion inset ----------------
ax = axs[5]
cap = d51["capacity"]
lnI = np.log(cap["doses"])
ax.errorbar(lnI, cap["N_means"], yerr=cap["N_sds"], fmt="o-", color=OI["blue"], ms=3.5,
            elinewidth=0.6, capsize=1.5, label="N(I0), MC ensemble")
cf = np.polyfit(lnI, cap["N_means"], 1)
ax.plot(lnI, np.polyval(cf, lnI), "--", color="0.4", lw=0.8, label=f"log law, slope = {cf[0]:.1f}")
ax.set_xlabel(r"ln $I_0$"); ax.set_ylabel("pulse count N")
ax.set_title("counting law; I(D;N) = 3.07 / 3.17 bits (96.8%)")
ax.legend(loc="upper left", fontsize=5.5)
axin = ax.inset_axes([0.52, 0.12, 0.45, 0.42])
g = d51["gain"]
gk = [r["sigma_k"] for r in g]
axin.plot(gk, [r["cv_A"] for r in g], "s-", color=OI["verm"], ms=2.5, label="CV(A)")
axin.plot(gk, [r["cv_T"] for r in g], "o-", color=OI["green"], ms=2.5, label="CV(T)")
axin.set_xlabel(r"$\sigma_k$", fontsize=5.5); axin.set_ylabel("CV", fontsize=5.5)
axin.tick_params(labelsize=5)
axin.set_title("gain-dispersion kill", fontsize=5.5)
axin.legend(fontsize=4.5, loc="upper left")
for spine in axin.spines.values():
    spine.set_linewidth(0.5)

for i, ax in enumerate(axs):
    panel_label(ax, "abcdef"[i])

fig.suptitle("EDFig. 4  p53 channel census: six-panel Fisher audit (code 51r3, M = 400 central, h = 2%)", fontsize=8)
fig.tight_layout(rect=[0, 0, 1, 0.96])
save(fig, os.path.join(OUT, "EDFig4_p53_channel_census_51r3"))
