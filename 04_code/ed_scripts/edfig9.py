import json, sys, os
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from edstyle import setup, panel_label, save, OI
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

setup()
ROOT = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911"
OUT = os.path.join(ROOT, r"05_主线纲领与设计\论文_细胞通讯审计_2026-09-23\figures_svg\ED")
P234 = os.path.join(ROOT, r"01_细胞线\结果\P2P3P4_盲裁决_本机原件")
P5 = os.path.join(ROOT, r"01_细胞线\结果\Keshelava2018_P5_47b勘误重裁")
P6 = os.path.join(ROOT, r"02_细胞线2\结果\P6_Wang2022_FCD归属")
P7 = os.path.join(ROOT, r"02_细胞线2\结果\P7_Moore2024_FCD-Weber")

VC = {"hit": OI["green"], "inter": OI["orange"], "fals": OI["verm"], "limit": OI["grey"]}

fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.6))
axs = axs.ravel()

# ---- a: clause tile ledger ----
ax = axs[0]
tiles = [  # (adjudication, [(clause, verdict)])
    ("P2", [("P2-1", "fals"), ("P2-2", "fals"), ("P2-3", "hit")]),
    ("P3", [("P3-1", "inter"), ("P3-2", "hit"), ("P3-3", "fals")]),
    ("P4", [("P4-1", "fals"), ("P4-2", "inter"), ("P4-3", "fals")]),
    ("P5", [("P5-1", "fals"), ("P5-2", "inter"), ("P5-3", "limit")]),
    ("P6", [("P6", "inter")]),
    ("P7", [("P7", "inter"), ("P7b", "inter")]),
]
y = 0
for adj, cl in tiles:
    for i, (c, v) in enumerate(cl):
        ax.add_patch(plt.Rectangle((i, y), 0.92, 0.85, color=VC[v], ec="white", lw=0.5))
        ax.text(i + 0.46, y + 0.42, c, ha="center", va="center", fontsize=5.5, color="white")
    ax.text(-0.15, y + 0.42, adj, ha="right", va="center", fontsize=7, fontweight="bold")
    y += 1
ax.set_xlim(-1.2, 3.1); ax.set_ylim(-0.4, y + 0.1)
ax.invert_yaxis(); ax.axis("off")
ax.set_title("frozen clause ledger: 14 clauses, 6 datasets")
handles = [Patch(color=VC["hit"], label="hit (2)"),
           Patch(color=VC["inter"], label="intermediate (5)"),
           Patch(color=VC["fals"], label="falsified (6)"),
           Patch(color=VC["limit"], label="dataset-level limitation (1)")]
ax.legend(handles=handles, loc="upper left", fontsize=5, ncol=2, bbox_to_anchor=(-0.28, -0.04))

# ---- b: P2 stratification units ----
ax = axs[1]
d44 = pd.read_csv(os.path.join(P234, "代码44b_CI单元表.csv"))
mol = d44[d44.channel == "mol_peak"].reset_index(drop=True)
tm = d44[d44.channel == "time_k3"].reset_index(drop=True)
xo = np.argsort(mol["delta"].values)
xx = np.arange(21)
ax.errorbar(xx, mol["delta"].values[xo],
            yerr=np.vstack([mol["delta"].values[xo] - mol["ci_lo"].values[xo],
                            mol["ci_hi"].values[xo] - mol["delta"].values[xo]]),
            fmt="o", ms=3, color=OI["blue"], elinewidth=0.5, capsize=0, label="mol channel (P2-1)")
ax.plot(xx, tm["delta"].values[xo], "s", ms=3, color=OI["green"], label="event time k=3 (P2-3)")
ax.axhline(0.03, color=OI["verm"], ls="--", lw=0.7)
ax.text(20.4, 0.033, "+0.03 hit line", fontsize=5, color=OI["verm"], ha="right")
ax.axhline(0, color="0.5", lw=0.5)
ax.axhspan(-0.01, 0.01, color=OI["green"], alpha=0.12)
ax.text(20.5, -0.0165, "P2-3 band |d| <= 0.01", fontsize=5, color=OI["green"], ha="right")
ax.set_xlabel("stratification unit (sorted, 7 promoters x 3 dose pairs)")
ax.set_ylabel("stratification penalty " + r"$\Delta$")
ax.set_title("P2 (Msn2): P2-1 0/7 pass, P2-3 6/7 pass")
ax.legend(fontsize=5.2, loc="upper left")
ax.set_xticks([0, 10, 20])

# ---- c: P3 AUCs ----
ax = axs[2]
d45 = pd.read_csv(os.path.join(P234, "代码45_单元表.csv"))
labels = ["dose 10>30", "dose 30>100", "mid>near", "far>mid", "dur 30>15", "dur 60>30"]
verds = ["inter", "inter", "hit", "hit", "fals", "fals"]
xx = np.arange(6)
for i, (_, r) in enumerate(d45.iterrows()):
    ax.errorbar(i, r["auc"], yerr=[[r["auc"] - r["ci_lo"]], [r["ci_hi"] - r["auc"]]],
                fmt="o", ms=4, color=VC[verds[i]], elinewidth=0.6, capsize=1.5)
ax.axhline(0.60, color="0.4", ls="--", lw=0.7); ax.axhline(0.55, color="0.4", ls=":", lw=0.7)
ax.text(5.4, 0.605, "0.60", fontsize=5, color="0.4", ha="right")
ax.text(5.4, 0.545, "0.55", fontsize=5, color="0.4", ha="right")
ax.axhline(0.5, color="0.6", lw=0.5)
ax.set_xticks(xx, labels, rotation=30, ha="right", fontsize=5.5)
ax.set_ylabel("AUC of free statistic")
ax.set_title("P3 (NF-kB gradient): 1 hit / 1 inter / 1 falsified")
ax.set_ylim(0.44, 1.0)

# ---- d: P4 AUCs + P4-3 inset ----
ax = axs[3]
d46 = pd.read_csv(os.path.join(P234, "代码46_单元表.csv"))
d412 = d46[d46.clause.isin(["P4-1", "P4-2"])].reset_index(drop=True)
labels = [f"{r['ligand'][:3]} {r['pair']}" for _, r in d412.iterrows()]
for i, (_, r) in enumerate(d412.iterrows()):
    ax.errorbar(i, r["auc"], yerr=[[r["auc"] - r["ci_lo"]], [r["ci_hi"] - r["auc"]]],
                fmt="o", ms=4, color=OI["verm"], elinewidth=0.6, capsize=1.5)
ax.axhline(0.60, color="0.4", ls="--", lw=0.7)
ax.text(0.1, 0.605, "0.60 line", fontsize=5, color="0.4")
ax.axhline(0.5, color="0.6", lw=0.5)
ax.set_xticks(range(8), labels, rotation=30, ha="right", fontsize=5)
ax.set_ylabel(r"$\tau(3)$ AUC")
ax.set_title("P4 (ERK-KTR): all 8 below 0.60")
ax.set_ylim(0.42, 0.72)
axin = ax.inset_axes([0.34, 0.63, 0.62, 0.33])
d43 = d46[d46.clause == "P4-3"].reset_index(drop=True)
cols = [OI["verm"] if abs(v) >= 0.03 else OI["green"] for v in d43["auc"]]
axin.bar(range(10), d43["auc"], color=cols, width=0.7)
axin.axhline(0.03, color="0.4", ls="--", lw=0.5); axin.axhline(-0.03, color="0.4", ls="--", lw=0.5)
axin.text(0.02, 0.90, "P4-3 size-gain stratification d", fontsize=4.5, transform=axin.transAxes, va="top")
axin.tick_params(labelsize=4.5)
axin.set_ylim(-0.16, 0.16)
for sp in axin.spines.values():
    sp.set_linewidth(0.4)

# ---- e: P5 AUCs ----
ax = axs[4]
d47 = pd.read_csv(os.path.join(P5, "代码47b_单元表.csv"))
statlab = {"ta": "event time", "cnt": "count", "pk": "paid peak"}
statcol = {"ta": OI["blue"], "cnt": OI["green"], "pk": OI["grey"]}
ax.axhspan(0.40, 0.60, color=OI["blue"], alpha=0.08)
ax.text(11.4, 0.42, "dead-zone band\n[0.40, 0.60)", fontsize=5, color=OI["blue"], ha="right", va="bottom")
for i, (_, r) in enumerate(d47.iterrows()):
    ax.errorbar(i, r["auc"], yerr=[[r["auc"] - r["lo"]], [r["hi"] - r["auc"]]],
                fmt="o", ms=4, color=statcol[r["stat"]], elinewidth=0.6, capsize=1.5)
ax.axhline(0.60, color="0.4", ls="--", lw=0.7)
abbr = {"ta": "tau", "cnt": "N", "pk": "pk"}
ax.set_xticks(range(len(d47)),
              [f"{abbr[r['stat']]} {r['pair'].replace(chr(0x2194), '-')}" for _, r in d47.iterrows()],
              rotation=45, ha="right", fontsize=4.5)
ax.set_ylabel("AUC")
ax.set_title("P5 (GPCR-Ca2+, 47b): dead-zone falsified")
ax.set_ylim(0.36, 0.85)

# ---- f: P6/P7 forest plot ----
ax = axs[5]
rows = [("P6 (NF-kB seq., 18 units)", -0.0014, -0.2518, 0.1533, "inter"),
        ("P7 first pass (25 units)", 0.1181, -0.2956, 0.4677, "inter"),
        ("P7b erratum (30 units)", -0.0909, -0.4446, 0.2766, "inter")]
for i, (lab, v, lo, hi, vd) in enumerate(rows):
    ax.errorbar(v, i, xerr=[[v - lo], [hi - v]], fmt="o", ms=5, color=VC[vd],
                elinewidth=0.8, capsize=2)
ax.axvline(0.15, color=OI["green"], ls="--", lw=0.7)
ax.axvline(-0.15, color=OI["verm"], ls="--", lw=0.7)
ax.axvline(0, color="0.5", lw=0.5)
ax.text(0.16, 2.3, "+0.15 hit line", fontsize=5, color=OI["green"])
ax.text(-0.16, 2.3, "-0.15 falsification line", fontsize=5, color=OI["verm"], ha="right")
ax.set_yticks(range(3), [r[0] for r in rows], fontsize=6)
ax.set_xlabel(r"$\Delta R^2$ = $R^2_{fold}$ - $R^2_{abs}$")
ax.set_title("P6 / P7 attribution: CIs cross zero")
ax.set_xlim(-0.55, 0.6); ax.set_ylim(-0.6, 2.7)
ax.invert_yaxis()

for i, ax in enumerate(axs):
    panel_label(ax, "abcdef"[i], x=-0.10)

fig.suptitle("EDFig. 9  Blind adjudications P2-P7: per-clause panels against frozen decision lines "
             "(seeds 20260814/20260815, bootstrap x2000)", fontsize=8)
fig.tight_layout(rect=[0, 0, 1, 0.955])
save(fig, os.path.join(OUT, "EDFig9_blind_adjudications_P2_P7"))
