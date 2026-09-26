# -*- coding: utf-8 -*-
"""
Fig. 6 | The pharmacological meta-audit.
Numbers transcribed from _archive_旧版本/Nature_SI_v04_2026-09-23.md S9, code 63/65/66 run outputs,
code80_arm_coordinate_residuals.csv, and atlas_audit CSVs (pair deviations
recomputed from atlas_audit_entries.csv, reproducing code 67: 17,987 pairs,
21.5% > 0.3 dex, 4.2% > 1 dex).
"""
import itertools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

plt.rcParams.update({
    "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "svg.fonttype": "none",
    "axes.unicode_minus": False,
})
BLUE, VERM, GREEN, ORANGE, GREY, SKY = ("#0072B2", "#D55E00", "#009E73",
                                        "#E69F00", "#7F7F7F", "#56B4E9")

ROOT = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911"
RES = ROOT + r"\04_细胞线4\_归档\结果"
OUT = ROOT + r"\05_主线纲领与设计\论文_细胞通讯审计_2026-09-23\figures_svg"

# --- load archived data ---
ent = pd.read_csv(RES + r"\atlas_audit\atlas_audit_entries.csv")
cols = ["doi", "receptor", "Measured process", "Pathway level", "Cell line",
        "Primary effector subtype"]
ent["grp"] = ent[cols].fillna("nan").astype(str).agg("|".join, axis=1)
pairs = []
for k, v in ent.groupby("grp"):
    d = v.drop_duplicates("ligand")
    if len(d) < 2:
        continue
    for i, j in itertools.combinations(d.index, 2):
        pairs.append(abs((d.loc[i, "tc"] - d.loc[j, "tc"])
                         - (d.loc[i, "logRA"] - d.loc[j, "logRA"])))
pairs = np.array(pairs)
assert len(pairs) == 17987, len(pairs)

pap = pd.read_csv(RES + r"\atlas_audit\atlas_audit_papers.csv")
r80 = pd.read_csv(RES + r"\code80_arm_coordinate_residuals.csv")

fig = plt.figure(figsize=(7.2, 5.6))
gs = gridspec.GridSpec(2, 2, height_ratios=[1.0, 1.12], hspace=0.62, wspace=0.32,
                       left=0.075, right=0.975, top=0.94, bottom=0.10)

def panel_label(ax, s, dx=-0.18, dy=1.06):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="top", ha="left")

# ---------------------------------------------------------------- panel a
# D2R fracture map (code 63 v1.1.0 run output; SI S9.2)
axa = fig.add_subplot(gs[0, 0])
ligs = ["Dopamine", "Ropinirole", "Aripiprazole", "Cariprazine",
        "Bifeprunox", "Pardoprunox", "S-3PPP"]
# (dCAMP, sdCAMP, dGaoB, sdGaoB); None = tau->inf, unresolvable
d2r = {
    "Dopamine":     (1.32, 0.45, 1.19, 0.15),
    "Ropinirole":   (None, None, None, None),
    "Aripiprazole": (-2.57, 0.10, -2.88, 0.16),
    "Cariprazine":  (-1.15, 0.11, -1.36, 0.15),
    "Bifeprunox":   (-3.54, 0.16, -3.46, 0.17),
    "Pardoprunox":  (0.18, 0.10, -0.51, 0.16),
    "S-3PPP":       (0.38, 0.10, 0.23, 0.14),
}
load_bearing = {"Aripiprazole", "Cariprazine", "Bifeprunox"}
x = np.arange(len(ligs))
for i, lig in enumerate(ligs):
    dc, sdc, dg, sdg = d2r[lig]
    if dc is None:
        axa.text(i, 0.55, "tau -> inf\n(unresolvable)", fontsize=5.2,
                 ha="center", va="center", color=GREY)
        continue
    col = VERM if lig == "Bifeprunox" else (ORANGE if lig in load_bearing else BLUE)
    axa.errorbar(i - 0.12, dc, yerr=sdc, fmt="o", ms=4.5, mfc="white",
                 color=col, mew=1.2, elinewidth=0.7, capsize=2, zorder=3)
    axa.errorbar(i + 0.12, dg, yerr=sdg, fmt="s", ms=4, mfc="white",
                 color=col, mew=1.2, elinewidth=0.7, capsize=2, zorder=3)
axa.axhline(0, color="black", lw=0.6)
axa.annotate("bifeprunox: -3.5 dex, 21-23 sigma", (4.12, -3.5),
             xytext=(1.1, -3.15), fontsize=6, color=VERM,
             arrowprops=dict(arrowstyle="-", lw=0.5, color=VERM))
axa.text(3.0, 1.55, "load-bearing for the study's\nsignificant-bias claims",
         fontsize=5.6, color=ORANGE, ha="center")
for i, lig in enumerate(ligs):
    if lig in load_bearing:
        axa.plot(i, 1.28, marker="v", ms=3.5, color=ORANGE, clip_on=False)
axa.text(0.02, -0.24, "post-repair bifeprunox residual 0.9-2.6 dex (unresolved)\n"
         "functional $pK_A$ conserved across arms (0.08-0.31 dex)",
         transform=axa.transAxes, fontsize=5.6, color="0.25", va="top")
axa.set_xticks(x)
axa.set_xticklabels([l if l != "Ropinirole" else "Ropinirole" for l in ligs],
                    fontsize=5.6, rotation=28, ha="right")
axa.set_xlim(-0.6, 6.6); axa.set_ylim(-4.0, 2.0)
axa.set_ylabel("$\\Delta$ = $pK_A^{func}$ - $pK_i$ (dex)")
axa.set_title("D2R fracture map (Klein Herenbrink 2016)", fontsize=7)
axa.errorbar([], [], yerr=1, fmt="o", ms=4.5, mfc="white", color=BLUE,
             label="cAMP arm")
axa.errorbar([], [], yerr=1, fmt="s", ms=4, mfc="white", color=BLUE,
             label="GαoB arm")
axa.legend(fontsize=5.6, loc="lower left", frameon=False, ncol=2,
           bbox_to_anchor=(0.0, 0.02))
panel_label(axa, "a")

# ---------------------------------------------------------------- panel b
# AT1R cross-arm fracture (code 65 run output; SI S9.3)
axb = fig.add_subplot(gs[0, 1])
# x positions: WT (0) and L112A (1); three readings per condition
conds = ["WT", "L112A"]
data = {
    "WT":    {"Gq": (6.11, 0.40), "bind": (7.50, 0.02), "arr": (8.16, 0.10),
              "d": -2.05, "sig": "5.0"},
    "L112A": {"Gq": (6.16, 0.30), "bind": (7.50, 0.02), "arr": (7.69, 0.08),
              "d": -1.53, "sig": "4.9"},
}
offs = {"Gq": -0.22, "bind": 0.0, "arr": 0.22}
style = {"Gq": ("o", BLUE, "Gq arm $pK_A$"),
         "bind": ("D", "black", "binding $pK_i$ = 7.50"),
         "arr": ("s", VERM, "β-arrestin arm $pK_A$")}
for j, cnd in enumerate(conds):
    for k in ["Gq", "bind", "arr"]:
        v, sd = data[cnd][k]
        m, c, _ = style[k]
        axb.errorbar(j + offs[k], v, yerr=sd, fmt=m, ms=5, mfc="white",
                     color=c, mew=1.2, elinewidth=0.7, capsize=2, zorder=3)
    axb.annotate("", xy=(j + 0.22, data[cnd]["arr"][0]),
                 xytext=(j - 0.22, data[cnd]["Gq"][0]),
                 arrowprops=dict(arrowstyle="<->", lw=0.6, color="0.4",
                                 shrinkA=3, shrinkB=3))
    axb.text(j + 0.34, (data[cnd]["arr"][0] + data[cnd]["Gq"][0]) / 2,
             f"Δ = {data[cnd]['d']:.2f} dex\n({data[cnd]['sig']}σ)",
             fontsize=6, color="0.25", ha="left", va="center")
for k, (m, c, lab) in style.items():
    axb.plot([], [], m, color=c, ms=5, mfc="white", mew=1.2, label=lab)
axb.legend(fontsize=5.2, loc="upper right", frameon=False, handlelength=1.2,
           borderaxespad=0.2)
axb.text(0.5, 5.62, "TRV055: Δ ≥ 1.6 dex (full-agonist lower bound, WT)",
         fontsize=5.6, color="0.25", ha="center")
axb.text(-0.45, 8.68, "binding affinity sits between the two arm readings",
         fontsize=5.6, color="0.25", ha="left", va="top")
axb.set_xticks([0, 1]); axb.set_xticklabels(["TRV026, WT", "TRV026, L112A"])
axb.set_xlim(-0.55, 1.75); axb.set_ylim(5.4, 8.8)
axb.set_ylabel("affinity (dex, $pK_A$ or $pK_i$)")
axb.set_title("AT1R cross-arm fracture (Wingler 2020)", fontsize=7)
panel_label(axb, "b")

# ---------------------------------------------------------------- panel c
# muOR GIRK-column fracture and code-80 repair (code 66 / code 80; SI S9.4, S9.8)
gsc = gs[1, 0].subgridspec(1, 2, width_ratios=[1, 1.15], wspace=0.55)
axc1 = fig.add_subplot(gsc[0])
# identity-1 GIRK column per ligand (code 66 run output)
girk = [("Morphine", 0.78, 4.5), ("Oliceridine", 0.69, 6.7),
        ("PZM21", 0.97, 8.4), ("SR-17018", 0.83, 3.9),
        ("Buprenorphine", 0.68, 5.1)]
names = [g[0] for g in girk]; vals = [g[1] for g in girk]; sigs = [g[2] for g in girk]
ypos = np.arange(len(girk))[::-1]
axc1.barh(ypos, vals, height=0.62, color=VERM, alpha=0.85, zorder=2)
for y, v, s, nm in zip(ypos, vals, sigs, names):
    axc1.text(0.02, y + 0.44, nm, fontsize=5.2, va="bottom", ha="left",
              color="0.15", zorder=4)
    axc1.text(v + 0.03, y, f"{v:.2f} ({s}σ)", fontsize=5.4, va="center",
              color="0.25")
axc1.axvline(0.3, color="black", lw=0.7, ls="--")
axc1.text(0.31, 5.0, "0.3 dex", fontsize=5.2, ha="left", va="bottom")
axc1.set_yticks(ypos); axc1.set_yticklabels([])
axc1.set_xlim(0, 1.35); axc1.set_ylim(-0.6, 5.5)
axc1.set_xlabel("|Δ| identity 1 (dex)", fontsize=6)
axc1.set_title("GIRK column fracture, per ligand\n(cross-assay map up to 12.4σ)",
               fontsize=6, pad=3)
panel_label(axc1, "c", dx=-0.18)

axc2 = fig.add_subplot(gsc[1])
# code-80 repair: |residual| before (M0) vs after (M1), archived CSV
m0 = r80["residual"].abs().to_numpy()
m1 = r80["residual_m1"].abs().to_numpy()
gmask = (r80["arm"] == "GIRK").to_numpy()
rng = np.random.default_rng(0)
for arr, xc, col in [(m0, 0, GREY), (m1, 1, GREY)]:
    axc2.scatter(np.full((~gmask).sum(), xc) + rng.uniform(-0.09, 0.09, (~gmask).sum()),
                 arr[~gmask], s=4, color=col, alpha=0.45, lw=0, zorder=2)
for arr, xc in [(m0, 0), (m1, 1)]:
    axc2.scatter(np.full(gmask.sum(), xc) + rng.uniform(-0.09, 0.09, gmask.sum()),
                 arr[gmask], s=9, color=VERM, alpha=0.9, lw=0, zorder=3)
for a, b in zip(m0[gmask], m1[gmask]):
    axc2.plot([0.06, 0.94], [a, b], color=VERM, lw=0.4, alpha=0.45, zorder=2)
axc2.set_xticks([0, 1]); axc2.set_xticklabels(["before (M0)", "after (M1)"],
                                              fontsize=6)
axc2.set_xlim(-0.5, 1.5); axc2.set_ylim(-0.05, 1.75)
axc2.set_ylabel("|residual| (dex)", fontsize=6)
axc2.set_title("code-80 repair: one coordinate per arm", fontsize=6, pad=3)
axc2.text(0.5, 0.965, "GIRK median 0.811 -> 0.107 dex\n"
          "whole table 0.140 -> 0.065 dex\nΔlogML = +128.2 (66 residuals)",
          fontsize=5.2, ha="center", va="top", color="0.25",
          transform=axc2.transAxes)
axc2.scatter([], [], s=9, color=VERM, lw=0, label="GIRK rows")
axc2.scatter([], [], s=5, color=GREY, lw=0, label="other rows")
axc2.legend(fontsize=5.2, loc="center left", frameon=False,
            bbox_to_anchor=(0.02, 0.62), handletextpad=0.2)

# ---------------------------------------------------------------- panel d
# Atlas-wide violations (code 67/71; SI S9.5)
gsd = gs[1, 1].subgridspec(1, 2, width_ratios=[1.25, 1], wspace=0.5)
axd1 = fig.add_subplot(gsd[0])
bins = np.concatenate([np.arange(0, 2.0, 0.1), np.arange(2.0, 5.2, 0.4)])
n, edges, patches = axd1.hist(pairs, bins=bins, color=BLUE, alpha=0.85, zorder=2)
for p, e in zip(patches, edges[:-1]):
    if e >= 1.0:
        p.set_facecolor(VERM)
    elif e >= 0.3:
        p.set_facecolor(ORANGE)
axd1.set_yscale("log")
axd1.axvline(0.3, color="black", lw=0.7, ls="--")
axd1.axvline(1.0, color="black", lw=0.7, ls="--")
axd1.text(0.36, 16000, "0.3 dex: 21.5% of pairs", fontsize=5.4, va="top", ha="left")
axd1.text(1.08, 3000, "1 dex: 4.2%", fontsize=5.4, va="top", ha="left")
axd1.set_xlim(0, 5.2); axd1.set_ylim(0.8, 30000)
axd1.set_xlabel("identity violation |Δ| (dex)", fontsize=6)
axd1.set_ylabel("ligand pairs (count)", fontsize=6)
axd1.set_title("Atlas-wide violations\n(17,987 pairs, 269 groups)", fontsize=6, pad=3)
panel_label(axd1, "d", dx=-0.34)

axd2 = fig.add_subplot(gsd[1])
pm = np.sort(pap["dev_median"].to_numpy())
yp = np.arange(len(pm))
cols2 = np.where(pm > 0.3, VERM, GREY)
axd2.scatter(pm, yp, s=5, c=cols2, lw=0, zorder=3)
axd2.axvline(0.3, color="black", lw=0.7, ls="--")
axd2.set_xlim(0, 2.4); axd2.set_ylim(-2, 68)
axd2.set_xlabel("per-paper median |Δ| (dex)", fontsize=6)
axd2.set_ylabel("papers (ranked)", fontsize=6)
axd2.set_title("per-paper medians\n15 of 66 papers > 0.3 dex", fontsize=6, pad=3)
axd2.text(0.97, 0.03, "pipeline closure:\n3,931/3,931 exact;\nfractured GIRK column\n"
          "silently absent\nfrom the Atlas", fontsize=5.0, color="0.25",
          va="bottom", ha="right", transform=axd2.transAxes)

fig.savefig(OUT + r"\Fig6_gpcr_metaaudit.svg")
fig.savefig(OUT + r"\Fig6_gpcr_metaaudit.png", dpi=200)
print("saved Fig6")
