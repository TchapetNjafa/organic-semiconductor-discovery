#!/usr/bin/env python3
"""
Shared publication style for Paper-2 figures (J. Chem. Inf. Model.).
Adapted from Paper1_Digital_Discovery/figstyle.py for visual series-consistency
across the two-paper set, with Elsevier (not RSC) column widths.
Import this FIRST in every figure script: import figstyle_p2 as fs
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl

# Elsevier column widths (inches). Single column 9.0 cm, double column 19.0 cm
# (standard elsarticle guide-for-authors figures).
SINGLE = 9.0 / 2.54
DOUBLE = 19.0 / 2.54

# Same fixed palette as Paper 1 for series consistency.
COLORS = {
    "primary": "#2b6cb0",   # blue
    "valid":   "#2f855a",   # green — stable donors
    "artifact": "#c53030",  # red — reactive structures
    "accent":  "#d69e2e",   # gold
    "grey":    "#4a5568",
    "light":   "#cbd5e0",
    "camb":    "#2f855a",
    "wb97":    "#d69e2e",
    "stable2": "#38a169",   # lighter green — added stable donor
    "stable3": "#276749",   # dark green — added stable donor
    "stable4": "#3182ce",   # mid blue — added stable donor
}

# Colored by stability class: the four stable donors in green/blue,
# the three reactive structures in red/gold/orange.
MOL_COLORS = {
    "1712": COLORS["valid"],
    "9168": COLORS["stable2"],
    "17574": COLORS["stable3"],
    "18506": COLORS["stable4"],
    "4550": COLORS["artifact"],
    "17851": COLORS["accent"],
    "20778": "#dd6b20",     # orange — third reactive
}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 8.5,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "axes.linewidth": 0.6,
    "lines.linewidth": 1.2,
    "patch.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.frameon": False,
    "legend.handlelength": 1.4,
    "legend.columnspacing": 1.0,
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})


def save(fig, stem):
    for ext in ("pdf", "png"):
        fig.savefig(f"figures/{stem}.{ext}")
    print(f"wrote figures/{stem}.pdf / .png")
