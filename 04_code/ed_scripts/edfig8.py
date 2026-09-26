import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from edstyle import setup, panel_label, schematic_tag, save, OI
import matplotlib.pyplot as plt

setup()
ROOT = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911"
OUT = os.path.join(ROOT, r"05_主线纲领与设计\论文_细胞通讯审计_2026-09-23\figures_svg\ED")

fig = plt.figure(figsize=(7.2, 4.8))
gs = fig.add_gridspec(2, 2)

# ---- a: duality instance table (SI S6.6) ----
ax = fig.add_subplot(gs[0, 0])
ax.axis("off")
rows = [
    ["G0 rescaling (5/5)", "multiplicative", "waveform corr / freq", "1.0 / 1.0  invariant"],
    ["G2 drift 0.3x-3x (7pt)", "monotone reparam", "event-time degrad.", "0.0000%"],
    ["G2 drift 0.3x-3x (7pt)", "monotone reparam", "ampl. mean / peak", "72.9% / 200%"],
    ["Type-I (SNIC)", "gain sweep", "freq / amplitude", "graded (sqrt) / satur."],
    ["Type-II (subcr. Hopf)", "gain sweep", "freq / amplitude", "ungraded / ungraded"],
]
tbl = ax.table(cellText=rows,
               colLabels=["instance", "transform", "statistic", "archived result"],
               loc="center", cellLoc="left")
tbl.auto_set_font_size(False)
tbl.set_fontsize(5.0)
tbl.scale(1.0, 1.35)
for (r, c), cell in tbl.get_celld().items():
    cell.set_linewidth(0.4)
    if r == 0:
        cell.set_text_props(fontweight="bold")
    if r in (2, 3) and c == 3:
        cell.set_text_props(color=OI["verm"] if r == 3 else OI["green"])
ax.set_title("Morris-Lecar duality instances (SI S6.6)", fontsize=7.5, pad=8)

# ---- b: gain-drift scan: event time vs amplitude degradation ----
ax = fig.add_subplot(gs[0, 1])
cats = ["event time", "amplitude (mean)", "amplitude (peak)"]
vals = [0.0, 72.9, 200.0]
bars = ax.bar(cats, vals, width=0.55,
              color=[OI["green"], OI["verm"], OI["verm"]])
for b, v in zip(bars, vals):
    ax.text(b.get_x() + b.get_width() / 2, v + 4, f"{v:.4f}%" if v < 1 else f"{v:.1f}%",
            ha="center", fontsize=5.8)
ax.set_ylabel("degradation under gain drift (%)")
ax.set_ylim(0, 235)
ax.set_title("G2 gain drift 0.3x-3x: event time survives, amplitude does not")
ax.tick_params(axis="x", labelsize=6)
schematic_tag(ax)

# ---- c: Type-I / SNIC: frequency graded, amplitude saturated ----
ax = fig.add_subplot(gs[1, 0])
d = np.linspace(0, 1.5, 200)
freq = np.sqrt(np.maximum(d, 0)) * 0.8
amp = np.minimum(1.0, 0.55 + 0.9 * d)
ax.plot(d, freq, "-", color=OI["blue"], lw=1.2, label="frequency (sqrt law, graded)")
ax.plot(d, amp, "--", color=OI["verm"], lw=1.2, label="amplitude (saturates)")
ax.set_xlabel("drive above SNIC (au)"); ax.set_ylabel("normalized value")
ax.set_title("Type-I (SNIC): counting channel open, amplitude shut")
ax.legend(fontsize=5.2, loc="lower left")
schematic_tag(ax)

# ---- d: Type-II / subcritical Hopf: both channels shut ----
ax = fig.add_subplot(gs[1, 1])
d2 = np.linspace(0, 1.5, 200)
amp2 = np.where(d2 < 0.5, 0.02, 1.0 - 0.258 * 0.3 * np.log10(d2 / 0.5 + 1))
freq2 = np.where(d2 < 0.5, 0.02, 0.85 + 0.03 * (d2 - 0.5))
ax.plot(d2, amp2, "--", color=OI["verm"], lw=1.2, label="amplitude: jumps at onset, then flat")
ax.plot(d2, freq2, "-", color=OI["blue"], lw=1.2, label="frequency: finite at onset, flat")
ax.axvline(0.5, color="0.4", ls=":", lw=0.7)
ax.text(0.53, 0.12, "subcritical onset:\nlarge limit cycle", fontsize=5, color="0.35")
ax.text(0.55, 0.70, "control fit: a = -0.258", fontsize=5.5, color=OI["verm"])
ax.set_xlabel("drive above onset (au)"); ax.set_ylabel("normalized value")
ax.set_title("Type-II (subcritical Hopf): neither channel graded")
ax.legend(fontsize=5.2, loc="upper left")
ax.set_ylim(0, 1.15)
schematic_tag(ax)

panel_label(fig.axes[0], "a", x=-0.02, y=1.02)
panel_label(fig.axes[1], "b")
panel_label(fig.axes[2], "c")
panel_label(fig.axes[3], "d")

fig.suptitle("EDFig. 8  Morris-Lecar duality: event timing is the invariant channel, amplitude is not\n"
             "(archived audit values; panels b-d are schematic redraws)", fontsize=8)
fig.tight_layout(rect=[0, 0, 1, 0.91])
save(fig, os.path.join(OUT, "EDFig8_MorrisLecar_duality"))
