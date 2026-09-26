# Extended Data schematics for
# "Scale-free statistics are the currency of cellular communication"
# EDFig1  : mathematical framework schematic (SI S1)
# EDFig13 : audit-trail workflow (SI S10)
# Style: Arial/Helvetica fallback DejaVu Sans, base 7pt, svg.fonttype='none',
# colorblind-safe muted palette (Okabe-Ito), vector schematic, no chartjunk.

import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import (FancyBboxPatch, FancyArrowPatch, Circle,
                                Ellipse, Polygon, Rectangle, Arc)

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
from daimon_runtime import setup_plot
setup_plot()

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "svg.fonttype": "none",
    "mathtext.default": "regular",
})

OUT = Path(__file__).resolve().parent

# Okabe-Ito muted palette
BLUE   = "#0072B2"
SKY    = "#56B4E9"
GREEN  = "#009E73"
ORANGE = "#E69F00"
VERM   = "#D55E00"
GRAY   = "#4D4D4D"
LGRAY  = "#999999"


def rbox(ax, x, y, w, h, fc, ec, lw=0.8, r=1.0, alpha=1.0, ls="-", z=2):
    p = FancyBboxPatch((x, y), w, h,
                       boxstyle=f"round,pad=0,rounding_size={r}",
                       fc=fc, ec=ec, lw=lw, alpha=alpha, ls=ls, zorder=z,
                       mutation_aspect=1)
    ax.add_patch(p)
    return p


def arrow(ax, x0, y0, x1, y1, color=GRAY, lw=1.1, style="-|>", ms=7, ls="-",
          z=3, alpha=1.0):
    a = FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style,
                        mutation_scale=ms, color=color, lw=lw, ls=ls,
                        shrinkA=0, shrinkB=0, zorder=z, alpha=alpha)
    ax.add_patch(a)
    return a


def check_mark(ax, x, y, s=2.2, color=GREEN, lw=1.4, z=4):
    ax.plot([x - s, x - s * 0.25, x + s], [y, y - s * 0.9, y + s * 0.9],
            color=color, lw=lw, solid_capstyle="round", zorder=z)


def cross_mark(ax, x, y, s=1.8, color=VERM, lw=1.4, z=4):
    ax.plot([x - s, x + s], [y - s, y + s], color=color, lw=lw,
            solid_capstyle="round", zorder=z)
    ax.plot([x - s, x + s], [y + s, y - s], color=color, lw=lw,
            solid_capstyle="round", zorder=z)


# ----------------------------------------------------------------------------
# EDFig 1 : mathematical framework schematic
# ----------------------------------------------------------------------------
fig = plt.figure(figsize=(7.2, 4.6))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 180)
ax.set_ylim(0, 112)
ax.axis("off")

# ---------------- panel a : layered cascade + degeneracy group -------------
ax.text(3, 108.7, "a", fontsize=8, fontweight="bold", va="center")
ax.text(7, 108.7, "Cascade of layers and the degeneracy group $G_i$",
        fontsize=7, fontweight="bold", va="center")
ax.text(7, 104.8, "inset: an orbit of $G_i$-equivalent parameter sets "
        "collapses to a point in output space",
        fontsize=4.8, style="italic", color=GRAY, va="center")

layers = [
    ("Layer 1: receptor", 86.0, 100.0),
    ("Layer 2: encoder",  64.0,  78.0),
    ("Layer 3: readout",  42.0,  56.0),
]
LX0, LX1 = 4.0, 58.0
for name, y0, y1 in layers:
    ymid = (y0 + y1) / 2
    rbox(ax, LX0, y0, LX1 - LX0, y1 - y0, "white", BLUE, lw=1.0, r=1.2)
    ax.text(LX0 + 2.2, ymid + 4.2, name, fontsize=6, fontweight="bold",
            va="center")
    ax.text(LX0 + 2.2, ymid + 0.2, r"parameters $\theta_i \in \Theta_i$",
            fontsize=4.8, va="center", color=GRAY)
    ax.text(LX0 + 2.2, ymid - 3.6, r"$Y_i \sim P(\cdot \mid Y_{i-1}, \theta_i)$",
            fontsize=4.8, va="center", color=GRAY)
    # degeneracy-orbit inset
    ecx, ecy = 43.0, ymid
    e = Ellipse((ecx, ecy), 13.0, 8.0, fc=SKY, ec=BLUE, lw=0.9, alpha=0.25,
                zorder=3)
    ax.add_patch(e)
    for ang in (0.4, 2.2, 4.2):
        ax.add_patch(Circle((ecx + 6.5 * np.cos(ang), ecy + 4.0 * np.sin(ang)),
                            0.9, fc=BLUE, ec="none", zorder=4))
    arrow(ax, 50.5, ymid, 54.6, ymid, color=GRAY, lw=0.9, ms=6)
    ax.add_patch(Circle((55.8, ymid), 1.1, fc=VERM, ec="none", zorder=4))

# input arrow into Layer 1
arrow(ax, 20, 103.0, 20, 100.6, color=GRAY, lw=1.1)
ax.text(22, 101.8, "input $S(t)$", fontsize=5, va="center", color=GRAY)
# inter-layer arrows
arrow(ax, 26, 86.0, 26, 78.6, color=GRAY, lw=1.1)
ax.text(28, 82.2, "$Y_1$", fontsize=5, va="center", color=GRAY)
arrow(ax, 26, 64.0, 26, 56.6, color=GRAY, lw=1.1)
ax.text(28, 60.2, "$Y_2$", fontsize=5, va="center", color=GRAY)
# output arrow
arrow(ax, 26, 42.0, 26, 37.2, color=GRAY, lw=1.1)
ax.text(28, 39.4, "output $Y_3$", fontsize=5, va="center", color=GRAY)

# connectors from layers into panel b
for _, y0, y1 in layers:
    arrow(ax, LX1, (y0 + y1) / 2, 65.0, (y0 + y1) / 2, color=LGRAY, lw=0.7,
          ms=5, z=1)

# ---------------- panel b : survival set -----------------------------------
ax.text(66, 108.7, "b", fontsize=8, fontweight="bold", va="center")
ax.text(70, 108.7, "Survival set $T_i$ (interface algebra)", fontsize=7,
        fontweight="bold", va="center")
ax.text(70, 104.8, "generic interface $y = kx + b$, gain $k$ not identifiable",
        fontsize=4.8, style="italic", color=GRAY, va="center")

rows = [
    ("difference ratios", True),
    ("sign, order (ordinals)", True),
    ("counts", True),
    ("event times", True),
    ("amplitude, absolute scale", False),
    ("fold change ($b \\neq 0$)", False),
]
ry = 99.0
for label, ok in rows:
    ax.text(66.5, ry + 2.6, label, fontsize=5.3, va="center")
    if ok:
        arrow(ax, 67, ry, 92, ry, color=GREEN, lw=1.4, ms=8)
        check_mark(ax, 96.5, ry + 0.4)
        ax.text(99.5, ry, "survives", fontsize=5.3, va="center", color=GREEN)
    else:
        arrow(ax, 67, ry, 92, ry, color=LGRAY, lw=1.0, ms=7)
        cross_mark(ax, 79.5, ry)
        ax.text(96.5, ry, "killed by $G_i$", fontsize=5.3, va="center",
                color=VERM)
    ry -= 8.6
ax.text(66.5, 46.5, "amplitude survives only if the gain is pinned from outside",
        fontsize=4.8, style="italic", color=GRAY, va="center")

# ---------------- panel c : three-leg criterion -----------------------------
ax.text(118, 108.7, "c", fontsize=8, fontweight="bold", va="center")
ax.text(122, 108.7, "Three-leg criterion", fontsize=7, fontweight="bold",
        va="center")

pillars = [
    (120.5, "Leg 1", ["algebraic", "survival"], ["fail: scale leak"], SKY, BLUE),
    (139.5, "Leg 2", ["symbolic", "reachability"],
     ["fail: invariant but", "uninformative"], "#7CCBA2", GREEN),
    (158.5, "Leg 3", ["temporal", "bandwidth"],
     ["fail: erased", "by filtering"], "#F5C983", "#B87A00"),
]
for x0, leg, lines, fails, fc, ec in pillars:
    cx = x0 + 8.5
    for j, f in enumerate(fails):
        ax.text(cx, 105.2 - 3.0 * j, f, fontsize=4.6, style="italic",
                color=VERM, ha="center", va="center")
    rbox(ax, x0, 74.5, 17, 24.0, fc, ec, lw=0.9, r=0.8, alpha=0.45)
    ax.text(cx, 94.0, leg, fontsize=6, fontweight="bold", ha="center",
            va="center")
    for j, ln in enumerate(lines):
        ax.text(cx, 89.6 - 3.4 * j, ln, fontsize=5.2, ha="center", va="center")
ax.add_patch(Rectangle((120, 68.5), 56, 5.0, fc=GRAY, ec="none", zorder=2))
ax.text(148, 71.0, "statistic $X$ crosses the interface iff all three hold",
        fontsize=5.2, color="white", ha="center", va="center", zorder=3)

# ---------------- panel c lower : information-loss identity -----------------
ax.text(118, 63.5, "Information-loss identity", fontsize=6.3,
        fontweight="bold", va="center")
ax.text(148, 59.0, r"$H(S) = I(S;Y_n) + L_1 + L_2 + L_3$", fontsize=7.5,
        ha="center", va="center")

segs = [("I(S;Y_n)", 14.7, SKY, False), ("L_1", 9.5, BLUE, False),
        ("L_2", 7.8, BLUE, False), ("L_3", 19.9, VERM, True)]
bx, by, bh = 122.0, 49.0, 5.0
for lab, w, col, hot in segs:
    ax.add_patch(Rectangle((bx, by), w, bh, fc=col, ec="white", lw=0.6,
                           alpha=0.9 if hot else 0.55, zorder=2))
    ax.text(bx + w / 2, by - 2.4, f"${lab}$" if lab != "I(S;Y_n)"
            else "$I(S;Y_n)$", fontsize=4.8, ha="center", va="center",
            color=VERM if hot else GRAY)
    bx += w
ax.text(148, 42.6, r"weakest layer $i^{\ast} = \arg\max_i L_i$ "
        "(highlighted) sets the end-to-end bound", fontsize=4.8, ha="center",
        va="center", color=VERM)

# ---------------- panel d : three reference classes -------------------------
ax.text(3, 35.8, "d", fontsize=8, fontweight="bold", va="center")
ax.text(7, 35.8, "Three reference classes: who sets and maintains the "
        "reference", fontsize=7, fontweight="bold", va="center")

tiles = [
    (2.0, "Class I: calibration-free", GREEN,
     "marker: no reference; nothing to maintain",
     "difference ratios, signs, counts, event times"),
    (61.0, "Class II: structural constants", BLUE,
     "marker: set by selection, kept by physics",
     "$K_d$, $K_m$ as dimensional converters; no cellular payment"),
    (120.0, "Class III: paid references", "#B87A00",
     "marker: active feedback on the reference itself",
     "set points, internal standards; paid continuously, drifts if unmaintained"),
]
for x0, head, col, marker, body in tiles:
    rbox(ax, x0, 4.5, 56, 27, col, col, lw=0.8, r=1.0, alpha=0.10)
    ax.add_patch(Rectangle((x0 + 0.8, 25.5), 54.4, 5.2, fc=col, ec="none",
                           alpha=0.30, zorder=2))
    ax.text(x0 + 28, 28.1, head, fontsize=6, fontweight="bold", ha="center",
            va="center", zorder=3)
    ax.text(x0 + 2.5, 20.5, marker, fontsize=5.2, style="italic", va="center")
    ax.text(x0 + 2.5, 15.8, body, fontsize=5.2, va="center", color=GRAY)
    extra = {
        2.0: "no maintenance signature expected",
        61.0: "no maintenance signature; outside class-III jurisdiction",
        120.0: "two of three markers + reference-value specificity",
    }[x0]
    ax.text(x0 + 2.5, 11.1, extra, fontsize=4.8, va="center", color=GRAY)

fig.savefig(OUT / "EDFig1_framework.svg")
fig.savefig(OUT / "EDFig1_framework.png", dpi=200)
plt.close(fig)

# ----------------------------------------------------------------------------
# EDFig 13 : audit-trail workflow
# ----------------------------------------------------------------------------
fig = plt.figure(figsize=(7.2, 2.75))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 184)
ax.set_ylim(14, 87.5)
ax.axis("off")

# rule banner
ax.add_patch(Rectangle((2, 78.5), 180, 7, fc=GRAY, ec="none", alpha=0.13))
ax.text(92, 82, "Rule: every quantitative claim traces to an archived script "
        "output, not to memory", fontsize=6.5, fontweight="bold",
        color=GRAY, ha="center", va="center")

# ---------------- artefact chain nodes --------------------------------------
NW, NH, NY = 26.0, 16.0, 50.0
node_x = [3.0 + i * (NW + 3.4) for i in range(6)]
titles = ["Pre-registration", "Decision table", "Scoring script",
          "Verdict card", "Errata", "Cloud re-execution"]
rules = [
    "hypotheses, clause list,\ndecision lines, tolerances,\nseeds fixed before data open",
    "thresholds frozen before\nscoring; never edited\nafterwards",
    "fixed seed\nSEED = 20260815\nbootstrap x2000",
    "frozen at first execution;\nverdicts never move",
    "corrections appended as\n44b / 47b / 49b;\nnever overwrite",
    "independent environment\nfrom raw public data;\nbit-identical check",
]


def icon_doc(cx, cy):
    ax.add_patch(Rectangle((cx - 3.0, cy - 2.6), 6.0, 5.6, fc="white",
                           ec=GRAY, lw=0.8, zorder=4))
    ax.add_patch(Polygon([(cx + 1.2, cy + 3.0), (cx + 3.0, cy + 3.0),
                          (cx + 3.0, cy + 1.2)], closed=True, fc="white",
                         ec=GRAY, lw=0.8, zorder=4))
    ax.plot([cx - 1.8, cx + 0.6], [cy + 0.6, cy + 0.6], color=LGRAY, lw=0.7,
            zorder=5)
    ax.plot([cx - 1.8, cx + 0.6], [cy - 0.7, cy - 0.7], color=LGRAY, lw=0.7,
            zorder=5)


def icon_table(cx, cy):
    ax.add_patch(Rectangle((cx - 3.0, cy - 2.6), 6.0, 5.6, fc="white",
                           ec=GRAY, lw=0.8, zorder=4))
    ax.plot([cx, cx], [cy - 2.6, cy + 3.0], color=GRAY, lw=0.6, zorder=5)
    for yy in (cy - 0.8, cy + 1.0):
        ax.plot([cx - 3.0, cx + 3.0], [yy, yy], color=GRAY, lw=0.6, zorder=5)


def icon_gear(cx, cy):
    for ang in np.linspace(0, 2 * np.pi, 8, endpoint=False):
        ax.plot([cx + 2.0 * np.cos(ang), cx + 3.1 * np.cos(ang)],
                [cy + 2.0 * np.sin(ang), cy + 3.1 * np.sin(ang)],
                color=GRAY, lw=1.0, zorder=4)
    ax.add_patch(Circle((cx, cy), 2.0, fc="white", ec=GRAY, lw=0.9, zorder=4))
    ax.add_patch(Circle((cx, cy), 0.8, fc=GRAY, ec="none", zorder=5))


def icon_card(cx, cy):
    icon_doc(cx, cy)
    check_mark(ax, cx + 0.2, cy + 0.6, s=1.6, color=GREEN, lw=1.5, z=6)


def icon_errata(cx, cy):
    icon_doc(cx, cy)
    ax.add_patch(Rectangle((cx - 3.0, cy - 4.3), 6.0, 1.7, fc=VERM, ec=VERM,
                           lw=0.6, alpha=0.30, zorder=5))
    ax.text(cx, cy - 3.45, "+", fontsize=5, fontweight="bold", color=VERM,
            ha="center", va="center", zorder=6)


def icon_cloud(cx, cy):
    ax.add_patch(Ellipse((cx, cy - 0.6), 6.4, 3.4, fc="white", ec=GRAY,
                         lw=0.8, zorder=4))
    ax.add_patch(Arc((cx - 1.5, cy + 0.6), 3.4, 3.0, theta1=40, theta2=180,
                     color=GRAY, lw=0.8, zorder=4))
    ax.add_patch(Arc((cx + 1.2, cy + 0.8), 3.8, 3.4, theta1=0, theta2=140,
                     color=GRAY, lw=0.8, zorder=4))
    check_mark(ax, cx + 0.1, cy + 0.3, s=1.2, color=GREEN, lw=1.2, z=6)


icons = [icon_doc, icon_table, icon_gear, icon_card, icon_errata, icon_cloud]

for x0, ttl, rule, ic in zip(node_x, titles, rules, icons):
    cx = x0 + NW / 2
    rbox(ax, x0, NY, NW, NH, "white", GRAY, lw=0.9, r=1.0)
    ic(cx, NY + NH - 5.2)
    ax.text(cx, NY + 3.0, ttl, fontsize=5.6, fontweight="bold", ha="center",
            va="center")
    for j, ln in enumerate(rule.split("\n")):
        ax.text(cx, NY - 2.2 - 2.9 * j, ln, fontsize=4.8, ha="center",
                va="center", color=GRAY)

for i in range(5):
    arrow(ax, node_x[i] + NW, NY + NH / 2, node_x[i + 1], NY + NH / 2,
          color=GRAY, lw=1.2, ms=8)

# ---------------- failure-inventory funnel + ledger -------------------------
fboxes = [
    ("falsified clauses", "P2-1, P2-2, P3-3, P4-1 ...", 34.5),
    ("rejected patches", "NF-kB repair codes 12, 15", 26.5),
    ("pipeline errors", "chemotaxis codes 55, 57, 59", 18.5),
]
for head, sub, y0 in fboxes:
    rbox(ax, 3.0, y0, 28.0, 6.0, VERM, VERM, lw=0.7, r=0.7, alpha=0.10)
    ax.text(5.0, y0 + 4.1, head, fontsize=5.0, fontweight="bold",
            color=VERM, va="center")
    ax.text(5.0, y0 + 1.4, sub, fontsize=4.3, color=GRAY, va="center")
    arrow(ax, 31.0, y0 + 3.0, 35.2, 26.5, color=VERM, lw=0.7, ms=5, alpha=0.8)

ax.add_patch(Polygon([(35.5, 36.5), (35.5, 16.5), (48.0, 24.2),
                      (48.0, 28.8)], closed=True, fc=ORANGE, ec="#B87A00",
                     lw=0.8, alpha=0.30, zorder=2))
arrow(ax, 48.0, 26.5, 59.6, 26.5, color=GRAY, lw=1.1, ms=7)

rbox(ax, 60.0, 20.0, 118.0, 12.0, GRAY, GRAY, lw=0.8, r=1.0, alpha=0.08)
ax.text(62.5, 28.2, "Double-recorded failure ledger", fontsize=6.0,
        fontweight="bold", color=GRAY, va="center")
ax.text(62.5, 23.6, "failures archived with the same discipline as hits; "
        "nothing deleted, corrections appended", fontsize=5.0, color=GRAY,
        va="center")

# dashed route: ledger feeds the same archive as the chain
ax.plot([178, 181.8, 181.8], [26, 26, NY + NH / 2], color=LGRAY, lw=0.8,
        ls=(0, (3, 2)), zorder=1)
arrow(ax, 181.8, NY + NH / 2, node_x[5] + NW + 0.4, NY + NH / 2, color=LGRAY,
      lw=0.8, ms=6, z=1)
ax.text(183.4, 42, "feeds the same archive", fontsize=4.5, color=GRAY,
        rotation=90, ha="center", va="center")

fig.savefig(OUT / "EDFig13_workflow.svg")
fig.savefig(OUT / "EDFig13_workflow.png", dpi=200)
plt.close(fig)

print("done")
