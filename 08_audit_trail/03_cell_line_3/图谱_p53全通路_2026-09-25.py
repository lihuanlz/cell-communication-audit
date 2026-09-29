# -*- coding: utf-8 -*-
"""p53 全通路上下游蛋白/酶链路图谱 v2（2026-09-25）。
按 L0-L4 级联层排布真实分子，颜色=单细胞活细胞时序数据可得性。
输出：03_细胞线3/结果/图谱_p53全通路上下游蛋白与酶_2026-09-25.{png,svg}
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patches as mpatches
from pathlib import Path
import sys

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
try:
    from daimon_runtime import setup_plot
    setup_plot()
except Exception:
    pass
plt.rcParams["svg.fonttype"] = "none"

ROOT = Path(__file__).resolve().parent.parent
OUT_PNG = ROOT / "结果" / "图谱_p53全通路上下游蛋白与酶_2026-09-25.png"
OUT_SVG = ROOT / "结果" / "图谱_p53全通路上下游蛋白与酶_2026-09-25.svg"

C_GREEN, E_GREEN = "#d9f0e3", "#0b6b3a"
C_YEL, E_YEL = "#fdf1dc", "#b07a1e"
C_RED, E_RED = "#fbe4e0", "#a03123"
GREEN, YEL, RED = "g", "y", "r"
FC = {GREEN: C_GREEN, YEL: C_YEL, RED: C_RED}
EC = {GREEN: E_GREEN, YEL: E_YEL, RED: E_RED}

# 六列：L0, L1, L1.5, L2, L3, L4
COLS = [0.005 + i * 0.166 for i in range(6)]
COLW = 0.152
BOX_H = 0.135
TOP = 0.86
BOT = 0.16

layers = [
    ("L0 damage", [
        ("DSB (35/Gy)", "lesion; Poisson\ncounting", GREEN),
        ("gammaH2AX foci", "amplifier / marker\n(H2AX S139-P)", YEL),
    ]),
    ("L1 sensing", [
        ("MRN complex\nMRE11-RAD50-NBS1", "DSB sensor,\nrecruits ATM", RED),
        ("Ku70/Ku80\n+ DNA-PKcs", "alt. sensor (NHEJ)", RED),
        ("ATM (kinase)", "apical kinase,\nS1981 auto-P", YEL),
        ("ATR + ATRIP\n(+ TOPBP1)", "ssDNA / replication\nstress branch", RED),
        ("53BP1 / MDC1\n/ BRCA1", "mediators", RED),
    ]),
    ("L1.5 transduction", [
        ("Chk2 (kinase)", "ATM->p53\n(p53 S15/20-P)", RED),
        ("Chk1 (kinase)", "ATR->p53", RED),
        ("p38 MAPK", "stress input to p53", RED),
        ("CDC25A/B/C\n(phosphatase)", "fast checkpoint\n(cell cycle brake)", RED),
    ]),
    ("L2 core oscillator", [
        ("p53", "transcription factor;\nstabilized by S15-P", GREEN),
        ("MDM2", "E3 ubiquitin ligase,\nnegative feedback", GREEN),
        ("MDM4", "p53-binding inhibitor", RED),
        ("Wip1 / PPM1D\n(PP2C phosphatase)", "master reset: de-P of ATM,\nChk1/2, p53-S15, MDM2,\ngammaH2AX", YEL),
    ]),
    ("L3 decoding (targets)", [
        ("p21 / CDKN1A", "arrest; leaky-integrator\ndecoder", GREEN),
        ("PUMA / BBC3", "apoptosis accumulator\n(high threshold)", YEL),
        ("BAX / NOXA", "apoptosis execution", RED),
        ("GADD45A /\n14-3-3sigma", "repair / G2 arrest", RED),
    ]),
    ("L4 fate", [
        ("arrest", "CDK-cyclin inhibition\n(via p21)", GREEN),
        ("apoptosis", "caspase-9/3, APAF1,\ncytochrome c", YEL),
        ("senescence", "permanent arrest", RED),
    ]),
]

fig, ax = plt.subplots(figsize=(20, 9.2))
ax.set_xlim(0, 1); ax.set_ylim(-0.30, 1.06); ax.axis("off")

def slot_y(i, n):
    if n == 1:
        return (TOP + BOT) / 2
    return TOP - i * (TOP - BOT) / (n - 1)

centers = {}
for ci, (name, prots) in enumerate(layers):
    x0 = COLS[ci]
    ax.text(x0 + COLW / 2, 1.005, name, ha="center", va="bottom", fontsize=11.5,
            fontweight="bold", color="#222222")
    n = len(prots)
    for i, (pname, role, lev) in enumerate(prots):
        y0 = slot_y(i, n)
        ax.add_patch(FancyBboxPatch((x0, y0 - BOX_H / 2), COLW, BOX_H,
                                    boxstyle="round,pad=0.005", fc=FC[lev],
                                    ec=EC[lev], lw=1.6))
        ax.text(x0 + COLW / 2, y0 + BOX_H * 0.20, pname, ha="center", va="center",
                fontsize=9.3, fontweight="bold", color=EC[lev])
        ax.text(x0 + COLW / 2, y0 - BOX_H * 0.26, role, ha="center", va="center",
                fontsize=7.3, color="#444444")
        centers[(ci, i)] = (x0 + COLW / 2, y0)

# 层间主干箭头（中线）
ymid = slot_y(0, 1)
for ci in range(5):
    ax.add_patch(FancyArrowPatch((COLS[ci] + COLW + 0.004, ymid),
                                 (COLS[ci + 1] - 0.004, ymid),
                                 arrowstyle="-|>", mutation_scale=20,
                                 color="#33527a", lw=2.0, zorder=1))

# 反馈箭头：MDM2 -> p53（L2 列内，负反馈）
x_p53, y_p53 = centers[(3, 0)]
x_mdm2, y_mdm2 = centers[(3, 1)]
ax.add_patch(FancyArrowPatch((x_mdm2 + COLW * 0.42, y_mdm2 + BOX_H / 2),
                             (x_p53 + COLW * 0.42, y_p53 - BOX_H / 2),
                             arrowstyle="-|>", mutation_scale=14, color=E_RED,
                             lw=1.5, linestyle="--",
                             connectionstyle="arc3,rad=-0.35"))
ax.text(x_p53 + COLW * 0.62, (y_p53 + y_mdm2) / 2, "MDM2 degrades p53\n(negative feedback)",
        fontsize=7.2, color=E_RED, ha="left", style="italic")
# Wip1 -> ATM 复位（长弧，底部走线）
x_wip1, y_wip1 = centers[(3, 3)]
x_atm, y_atm = centers[(1, 2)]
ax.add_patch(FancyArrowPatch((x_wip1, y_wip1 - BOX_H / 2 - 0.005),
                             (x_atm, y_atm - BOX_H / 2 - 0.005),
                             arrowstyle="-|>", mutation_scale=14, color="#7a4a9e",
                             lw=1.5, linestyle="--",
                             connectionstyle="arc3,rad=0.18"))
ax.text((x_wip1 + x_atm) / 2, 0.045, "Wip1 resets the cascade (dephosphorylates ATM, Chk1/2, p53, MDM2, gammaH2AX)",
        fontsize=7.6, color="#7a4a9e", ha="center", style="italic")

# 审计判词横带
verdicts = [
    (0, "L0: Poisson floor\n0.28 bits (codes 81/83)", E_GREEN),
    (1, "L1: weakest layer 0.87 bits;\n98% mapping choice, not physics;\nno public ATM single-cell\nseries (registered gap)", E_RED),
    (2, "L1.5: folded into L1/L2\nin our cascade model;\nChk2 = ATM->p53 transducer", E_YEL),
    (3, "L2: near-lossless 0.012 bits;\nT=5.5 h pinned; Wip1 axis =\nstiffest Fisher direction\n(codes 81/82)", E_GREEN),
    (4, "L3: integrator decoders count,\ncannot read period\n(code 83)", E_GREEN),
    (5, "L4: fate threshold on\nPUMA/p21 ratio;\nend-to-end 1.33 bits", E_GREEN),
]
for ci, txt, col in verdicts:
    ax.text(COLS[ci] + COLW / 2, -0.045, txt, ha="center", va="top", fontsize=7.4,
            color=col,
            bbox=dict(boxstyle="round,pad=0.32", fc="white", ec=col, lw=0.9))

handles = [mpatches.Patch(fc=C_GREEN, ec=E_GREEN, label="public single-cell live time series available"),
           mpatches.Patch(fc=C_YEL, ec=E_YEL, label="reporter / method exists; public time series scarce"),
           mpatches.Patch(fc=C_RED, ec=E_RED, label="population-level or endpoint data only")]
ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.005, -0.30),
          fontsize=9.5, framealpha=0.95)
ax.set_title("p53 pathway: full protein/enzyme chain mapped onto the five-layer cascade "
             "(color = single-cell data availability)", fontsize=13.5, pad=26)

fig.savefig(OUT_PNG, dpi=200, bbox_inches="tight")
fig.savefig(OUT_SVG, bbox_inches="tight")
print("saved:", OUT_PNG)
