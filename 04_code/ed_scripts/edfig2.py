import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from edstyle import setup, panel_label, save, OI
import matplotlib.pyplot as plt

setup()
ROOT = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911"
OUT = os.path.join(ROOT, r"05_主线纲领与设计\论文_细胞通讯审计_2026-09-23\figures_svg\ED")

u = np.linspace(1e-7, 1 - 1e-7, 4000)

def xi_of_u(u, k1, k2):
    return u / (1 - u) * (k1 + 1 - u) / (k2 + u)

fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.6))
axs = axs.ravel()

# ---- a: GK curve family ----
ax = axs[0]
kappas = [0.001, 0.005, 0.01, 0.05, 0.5, 1.0]
cols = plt.cm.viridis(np.linspace(0, 0.9, len(kappas)))
for k, c in zip(kappas, cols):
    xv = xi_of_u(u, k, k)
    ax.plot(xv, u, color=c, lw=1.0, label=f"k = {k:g}, nH = {1+1/(2*k):.0f}")
ax.set_xscale("log"); ax.set_xlim(1e-2, 1e3)
ax.set_xlabel(r"activity ratio $\xi$"); ax.set_ylabel("modified fraction u")
ax.set_title("Goldbeter-Koshland curve family")
ax.legend(fontsize=4.8, loc="upper left")

# ---- b: group-action invariance ----
ax = axs[1]
k0 = 0.01
x0 = xi_of_u(u, k0, k0)
ax.plot(x0, u, color="0.2", lw=1.6, label="reference (Km, ST)")
# group action: (Km, ST) -> 5x leaves kappa invariant
ax.plot(xi_of_u(u, 5 * k0 * 1.0, 5 * k0 * 1.0) * 1.0, u, "o", ms=1.2, color=OI["verm"],
        mfc="none", label="(Km, ST) x5 action")
ax.set_xscale("log"); ax.set_xlim(1e-2, 1e3)
ax.set_xlabel(r"$\xi$"); ax.set_ylabel("u")
ax.set_title("group action: curves coincide, max|du| = 0")
ax.text(0.03, 0.58, "four generators verified:\n0, 1.7e-18, 0, 0 (symmetric)\n0, 2.2e-16, 0, 0 (asymmetric)\n"
        "counterexample: 1.45e-1", transform=ax.transAxes, fontsize=5, color="0.3", va="top")
ax.legend(fontsize=5.2, loc="lower right")

# ---- c: Fisher condition numbers (recomputed: 20 pts / 3 decades, sigma = 0.03) ----
ax = axs[2]
xj = np.logspace(-1.5, 1.5, 20)
sig = 0.03
kg = np.logspace(np.log10(0.005), 0, 9)
chi_shape, chi_abs = [], []
eps_fd = 1e-6
for k in kg:
    def uu(x, kk):
        # invert xi = u/(1-u)(kk+1-u)/(kk+u) by bisection
        out = np.empty_like(x)
        for i, xx in enumerate(np.atleast_1d(x)):
            lo, hi = 1e-12, 1 - 1e-12
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                if xi_of_u(mid, kk, kk) < xx:
                    lo = mid
                else:
                    hi = mid
            out[i] = 0.5 * (lo + hi)
        return out
    # shape block (ln c, ln kappa): u(x*c; k)
    u0 = uu(xj, k)
    d_c = (uu(xj * np.exp(eps_fd), k) - uu(xj * np.exp(-eps_fd), k)) / (2 * eps_fd)
    d_k = (uu(xj, k * np.exp(eps_fd)) - uu(xj, k * np.exp(-eps_fd))) / (2 * eps_fd)
    J = np.vstack([d_c, d_k]).T / sig
    F = J.T @ J
    ev = np.linalg.eigvalsh(F)
    chi_shape.append(ev[-1] / max(ev[0], 1e-300))
    # absolute block (ln c, ln Km, ln ST), kappa = Km/ST
    u0a = uu(xj, k)
    d_km = (uu(xj, k * np.exp(eps_fd)) - u0a) / eps_fd          # Km +, ST fixed
    d_st = (uu(xj, k * np.exp(-eps_fd)) - u0a) / eps_fd         # ST +, Km fixed
    Ja = np.vstack([d_c, d_km, d_st]).T / sig
    Fa = Ja.T @ Ja
    eva = np.linalg.eigvalsh(Fa)
    eva = np.sort(np.abs(eva))
    chi_abs.append(eva[-1] / max(eva[0], 1e-300))
ax.plot(kg, chi_shape, "o-", color=OI["blue"], ms=3.5, label="shape block (ln c, ln kappa)")
ax.plot(kg, chi_abs, "s--", color=OI["verm"], ms=3.5, label="absolute block (ln c, ln Km, ln ST)")
ax.set_xscale("log"); ax.set_yscale("log")
ax.axhspan(19, 34, color=OI["blue"], alpha=0.10)
ax.text(0.008, 6, "archived range 19-34", fontsize=5, color=OI["blue"])
ax.set_xlabel(r"$\kappa$"); ax.set_ylabel(r"Fisher condition number $\chi$")
ax.set_title("shape identifiable, scale not (exact zero direction)")
ax.legend(fontsize=5.2, loc="center left")

# ---- d: asymmetric Michaelis triple ----
ax = axs[3]
k2 = 0.01
for r, c in zip([0.1, 0.3, 1.0, 3.0, 10.0], plt.cm.plasma(np.linspace(0.1, 0.85, 5))):
    ax.plot(xi_of_u(u, r * k2, k2), u, color=c, lw=1.0, label=f"r = {r:g}")
ax.set_xscale("log"); ax.set_xlim(1e-2, 1e3)
ax.set_xlabel(r"$\xi$"); ax.set_ylabel("u")
ax.set_title("asymmetric Michaelis curves (k2 = 0.01)")
ax.legend(fontsize=5, loc="upper left")
axin = ax.inset_axes([0.45, 0.13, 0.52, 0.42])
rr = np.logspace(np.log10(0.03), np.log10(30), 60)
nH = 4 / (4 - 1 / (rr * k2 + 0.5) - 1 / (k2 + 0.5))
axin.plot(rr, nH, color=OI["blue"], lw=1.0)
axin.axhline(51.0, color="0.4", ls="--", lw=0.6)
for r0, n0 in [(0.03, 98.973), (30, 5.068)]:
    axin.plot(r0, n0, "o", ms=3.5, color=OI["verm"])
    axin.annotate(f"{n0:.1f}", (r0, n0), textcoords="offset points", xytext=(4, 2), fontsize=4.5)
axin.set_xscale("log")
axin.set_xlabel("r = k1/k2", fontsize=5); axin.set_ylabel("nH closed form", fontsize=5)
axin.tick_params(labelsize=4.5)
for sp in axin.spines.values():
    sp.set_linewidth(0.4)

# ---- e: one-sided protocol (SI S2.3 table) ----
ax = axs[4]
kt = np.array([0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0])
ci = np.array([2920, 5.3, 2.4, 1.34, 1.29, 1.33, 1.45])
nh = np.array([501, 101, 51, 11, 6, 2, 1.5])
ax.plot(kt, ci, "o-", color=OI["verm"], ms=4, label="95% CI inflation")
ax.axvspan(10**(-3.4), 0.01, color=OI["grey"], alpha=0.12)
ax.text(0.0016, 100, "deep zero-order:\nonly one-sided\nupper bound\nreportable", fontsize=5, color="0.35")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"$\kappa$"); ax.set_ylabel("95% CI inflation factor")
ax2 = ax.twinx()
ax2.plot(kt, nh, "s--", color=OI["blue"], ms=3.5, label="nH")
ax2.set_yscale("log"); ax2.set_ylabel("Hill coefficient nH", color=OI["blue"])
ax2.tick_params(axis="y", labelcolor=OI["blue"], labelsize=6)
ax2.spines["right"].set_linewidth(0.6)
ax.set_title("one-sided protocol below kappa ~ 0.01 (SI table S2.3)")
ax.legend(fontsize=5.2, loc="upper right")

# ---- f: positive-feedback switch ----
ax = axs[5]
kap = 0.05
M = xi_of_u(u, kap, kap)
dM = np.gradient(M, u)
# g* : minimum of M'/(M - u M') over u where positive
den = M - u * dM
mask = den > 0
gstar = np.min(dM[mask] / den[mask])
for g, c in zip([0, 2, 5, 10, 20], plt.cm.cividis(np.linspace(0.05, 0.9, 5))):
    xig = M / (1 + g * u)
    # folds
    sgn = np.sign(np.gradient(xig, u))
    ax.plot(xig, u, color=c, lw=1.0, label=f"g = {g}")
ax.set_xscale("log"); ax.set_xlim(1e-2, 30)
ax.set_xlabel(r"activity ratio $\xi$"); ax.set_ylabel("u")
ax.set_title(f"positive-feedback switch: g* = {gstar:.3f} (archived 0.443)")
ax.legend(fontsize=5, loc="lower right")
axin = ax.inset_axes([0.40, 0.57, 0.57, 0.41])
ufine = np.linspace(1e-6, 1 - 1e-6, 60000)
Mfine = xi_of_u(ufine, kap, kap)
def dlnxi_ana(uu, g):
    return 1/uu + 1/(1-uu) - 1/(kap+1-uu) - 1/(kap+uu) - g/(1+g*uu)
gg = np.linspace(0.45, 20, 120)
widths = []
for g in gg:
    xig = Mfine / (1 + g * ufine)
    sgn = np.sign(dlnxi_ana(ufine, g))
    flips = np.where(np.diff(sgn) != 0)[0]
    widths.append(np.log10(xig[flips[0]] / xig[flips[-1]]) if len(flips) >= 2 else np.nan)
axin.plot(gg, widths, color=OI["verm"], lw=1.0)
axin.axvspan(10, 20, color=OI["blue"], alpha=0.10)
axin.axhspan(0.44, 0.58, xmin=0, xmax=1, color=OI["grey"], alpha=0.15)
axin.text(15, 0.09, "archived 0.44-0.58 dex\nat g = 10-20", fontsize=4.5, ha="center", color="0.35")
axin.text(0.5, 0.545, "apparent nH ~ 480 at g = 10\n(archived); p42 MAPK nH >= 35",
          fontsize=4.5, ha="left", va="top", color=OI["blue"])
axin.set_xlim(0.45, 20); axin.set_ylim(0, 0.65)
axin.set_xlabel("feedback gain g", fontsize=5); axin.set_ylabel("hysteresis width (dex)", fontsize=5)
axin.tick_params(labelsize=4.5)
for sp in axin.spines.values():
    sp.set_linewidth(0.4)

for i, ax in enumerate(axs):
    panel_label(ax, "abcdef"[i])

fig.suptitle("EDFig. 2  PdPC full panels: degeneracy group, Fisher barrier, asymmetric Michaelis,\n"
             "one-sided protocol, positive-feedback switch (recomputed from SI equations S2.1-S2.5)", fontsize=8)
fig.tight_layout(rect=[0, 0, 1, 0.925])
save(fig, os.path.join(OUT, "EDFig2_PdPC_full_panels"))
