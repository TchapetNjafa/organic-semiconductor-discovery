#!/usr/bin/env python3
"""
Shared publication style for all Paper-1 figures (RSC Digital Discovery).
Import this FIRST in every figure script:  import figstyle as fs
It sets global matplotlib rcParams for a single consistent look and exports
standard column widths, a fixed colour palette, and a save helper.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl

# RSC column widths (inches). Single column 8.3 cm, double column 17.1 cm.
SINGLE = 8.3 / 2.54
DOUBLE = 17.1 / 2.54

# Fixed palette used across every figure.
COLORS = {
    "primary": "#2b6cb0",   # blue  — B3LYP / primary series
    "valid":   "#2f855a",   # green — chemically valid / stable
    "artifact":"#c53030",   # red   — artifacts / unphysical / reactive
    "accent":  "#d69e2e",   # gold  — third series
    "grey":    "#4a5568",   # dark grey — neutral text / controls
    "light":   "#cbd5e0",   # light grey — background point cloud
    "camb":    "#2f855a",
    "wb97":    "#d69e2e",
}

mpl.rcParams.update({
    # fonts — one sans-serif family everywhere; editable text in the PDF
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 8.5,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7,
    "pdf.fonttype": 42,          # TrueType — text stays selectable/editable
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    # lines / axes
    "axes.linewidth": 0.6,
    "lines.linewidth": 1.2,
    "patch.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    # ticks
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    # legend
    "legend.frameon": False,
    "legend.handlelength": 1.4,
    "legend.columnspacing": 1.0,
    # output
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})


def save(fig, stem):
    """Save a figure as vector PDF + 600-dpi PNG under figures/."""
    for ext in ("pdf", "png"):
        fig.savefig(f"figures/{stem}.{ext}")
    print(f"wrote figures/{stem}.pdf / .png")
