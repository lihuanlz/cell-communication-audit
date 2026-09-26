"""Shared style for Extended Data figures (SVG, fonttype=none)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Okabe-Ito colorblind-safe palette
OI = {
    "black": "#000000", "orange": "#E69F00", "sky": "#56B4E9",
    "green": "#009E73", "yellow": "#F0E442", "blue": "#0072B2",
    "verm": "#D55E00", "pink": "#CC79A7", "grey": "#7F7F7F",
}

def setup():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.linewidth": 0.6,
        "axes.labelsize": 7,
        "axes.titlesize": 7.5,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "legend.fontsize": 6,
        "legend.frameon": False,
        "svg.fonttype": "none",
        "figure.dpi": 200,
        "savefig.dpi": 200,
        "lines.linewidth": 1.0,
        "lines.markersize": 3.5,
    })

def panel_label(ax, letter, x=-0.14, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=9, fontweight="bold", va="top", ha="left")

def schematic_tag(ax):
    ax.text(0.99, 0.01, "schematic redraw of archived result",
            transform=ax.transAxes, fontsize=5, color="0.55",
            ha="right", va="bottom", style="italic")

def save(fig, outbase):
    fig.savefig(outbase + ".svg", bbox_inches="tight")
    fig.savefig(outbase + ".png", bbox_inches="tight", dpi=200)
    plt.close(fig)
    print("wrote", outbase + ".svg / .png")
