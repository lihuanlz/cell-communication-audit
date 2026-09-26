import json, sys, os
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from edstyle import setup, panel_label, save, OI
import matplotlib.pyplot as plt

setup()
ROOT = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911"
AT = os.path.join(ROOT, r"04_细胞线4\_归档\结果\atlas_audit")
OUT = os.path.join(ROOT, r"05_主线纲领与设计\论文_细胞通讯审计_2026-09-23\figures_svg\ED")

entries = pd.read_csv(os.path.join(AT, "atlas_audit_entries.csv"))
groups = pd.read_csv(os.path.join(AT, "atlas_audit_groups.csv"))
papers = pd.read_csv(os.path.join(AT, "atlas_audit_papers.csv"))
repair = pd.read_csv(os.path.join(AT, "atlas_repair_projection.csv"))
frac = pd.read_csv(os.path.join(AT, "atlas_fracture_probability.csv"))
summary = json.load(open(os.path.join(AT, "atlas_audit_summary.json")))

fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.8))
axs = axs.ravel()

# ---- a: entry-level deviation distribution ----
ax = axs[0]
dev = entries["dev_entry"].abs().dropna().values
bins = np.logspace(-3, 1.2, 60)
ax.hist(dev, bins=bins, color=OI["blue"], alpha=0.85)
ax.axvline(0.3, color=OI["verm"], ls="--", lw=0.8)
ax.axvline(1.0, color=OI["verm"], ls=":", lw=0.8)
ax.text(0.31, ax.get_ylim()[1]*0.8, "0.3 dex", fontsize=5, color=OI["verm"])
ax.text(1.02, ax.get_ylim()[1]*0.6, "1 dex", fontsize=5, color=OI["verm"])
ax.set_xscale("log")
ax.set_xlabel(r"|$\Delta$log($\tau$/K$_A$) - $\Delta$log(E$_{max}$/EC$_{50}$)|  (dex, per entry)")
ax.set_ylabel("entries")
ax.set_title("violation distribution (entries; pairs 21.5% > 0.3 dex)")
q = summary["pair_dev_quantiles_dex"]
ax.text(0.03, 0.95, f"pair quantiles: P50 {q['P50']:.2f}, P90 {q['P90']:.2f},\nP99 {q['P99']:.2f} dex (17987 pairs)",
        transform=ax.transAxes, fontsize=5, va="top", color="0.35")

# ---- b: per-paper medians ----
ax = axs[1]
pm = papers["dev_median"].sort_values().values
ax.plot(np.arange(1, len(pm) + 1), pm, "o", ms=2.5, color=OI["blue"])
ax.axhline(0.3, color=OI["verm"], ls="--", lw=0.8)
ax.set_yscale("log")
ax.set_xlabel("paper (rank, 66 multi-ligand papers)")
ax.set_ylabel("median |deviation| (dex)")
ax.set_title("per-paper medians: 15/66 (23%) > 0.3 dex")
n_exc = (pm > 0.3).sum()
ax.text(2, 0.35, f"{n_exc} papers above line", fontsize=5, color=OI["verm"])

# ---- c: repair-cost distribution (code 77) ----
ax = axs[2]
cost = repair["cost_rms_dex"].dropna().values
bins = np.logspace(-3, 1.2, 50)
ax.hist(cost, bins=bins, color=OI["green"], alpha=0.85)
ax.axvline(0.072, color="0.3", ls="-", lw=0.8)
ax.text(0.075, ax.get_ylim()[1]*0.85, "median 0.072", fontsize=5, color="0.3")
ax.axvline(0.3, color=OI["verm"], ls="--", lw=0.8)
ax.text(0.28, ax.get_ylim()[1]*0.7, "24.2% of groups > 0.3 ", fontsize=5, color=OI["verm"], ha="right")
ax.set_xscale("log")
ax.set_xlabel("group repair cost, RMS (dex)")
ax.set_ylabel("assay groups")
ax.set_title("constraint-projection repair cost (269 groups)")

# ---- d: repair move vs fracture probability ----
ax = axs[3]
ax.plot(np.abs(frac["repair_move"]), frac["p_broken"], ".", ms=1.5, color=OI["pink"], alpha=0.4)
ax.set_xscale("log")
ax.set_xlabel("|repair move| per entry (dex)")
ax.set_ylabel("fracture probability")
ax.set_title("repair cost vs fracture prob. (Spearman 0.962)")
ax.text(0.03, 0.9, "15.0% of entries p > 0.9\nvs 14.5% with |dev| > 0.3 dex\nEM mixture: pi_h = 0.745, "
        "sigma_h = 0.089, sigma_b = 0.683 dex",
        transform=ax.transAxes, fontsize=5, va="top", color="0.35")

# ---- e: threshold robustness (code 72 M5) ----
ax = axs[4]
thr = [0.2, 0.3, 0.5, 0.8, 1.0]
grp = [39.03, 26.02, 14.87, 7.81, 5.20]
pap = [37.88, 22.73, 13.64, 6.06, 3.03]
ax.plot(thr, grp, "o-", color=OI["blue"], ms=4, label="groups (269)")
ax.plot(thr, pap, "s--", color=OI["verm"], ms=4, label="papers (66)")
for x, y in zip(thr, grp):
    ax.annotate(f"{y:.1f}", (x, y), textcoords="offset points", xytext=(0, 4), fontsize=5,
                color=OI["blue"], ha="center")
ax.set_xlabel("deviation threshold (dex)"); ax.set_ylabel("exceedance rate (%)")
ax.set_title("threshold robustness (code 72 M5)")
ax.legend(fontsize=5.5)

# ---- f: robustness triple (code 73) ----
ax = axs[5]
cats = ["neg. control\n(Hb)", "aggregation\nrho 0-0.5", "conservative\nSIG=0.10", "jackknife"]
vals = [0, 55.6, 28.2, 65.7]
cols = [OI["green"], OI["blue"], OI["sky"], OI["orange"]]
bars = ax.bar(range(4), vals, color=cols, width=0.62)
notes = ["0/5 alarms\np50 spread\n0.0065 dex", "134/241\ngroups", "68/241\ngroups", "65.7% flagged\nGIRK jack. = 1.0"]
for b, v, n_ in zip(bars, vals, notes):
    ax.text(b.get_x() + b.get_width()/2, v + 1.5, n_, ha="center", fontsize=4.8)
ax.set_xticks(range(4), cats, fontsize=5.5)
ax.set_ylabel("%")
ax.set_title("robustness triple (code 73)")
ax.set_ylim(0, 85)

for i, ax in enumerate(axs):
    panel_label(ax, "abcdef"[i], x=-0.10)

fig.suptitle("EDFig. 12  Biased Signaling Atlas meta-audit: violation distributions, per-paper medians, "
             "constraint-projection repair (codes 67, 71, 72, 73, 77)", fontsize=8)
fig.tight_layout(rect=[0, 0, 1, 0.955])
save(fig, os.path.join(OUT, "EDFig12_atlas_violations_repair"))
