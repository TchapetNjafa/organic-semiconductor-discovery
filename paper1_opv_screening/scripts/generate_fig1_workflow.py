#!/usr/bin/env python3
"""
Figure 1: screening + chemical-validity workflow (matplotlib flowchart).
Dependency-free. Emits figures/figure1_workflow.{pdf,png}.
Counts are the verified pipeline numbers (docs/REFRAME_PLAN_2026-07-06.md).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import figstyle as fs
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8,
})

# ---- palette -----------------------------------------------------------
# Each stage type gets a saturated "ink" colour (edges, arrows, headers)
# and a very light tint of the same hue for the fill, so boxes read as
# a coherent family rather than flat white cards.
INK = {
    "grey":  "#475569",   # neutral pipeline steps
    "blue":  "#1d4ed8",   # computational steps
    "red":   "#b91c1c",   # raw / unfiltered candidates
    "green": "#047857",   # validity gate / final result
    "gold":  "#b45309",   # final, single-molecule highlight
}
FILL = {
    "grey":  "#f1f5f9",
    "blue":  "#eff6ff",
    "red":   "#fef2f2",
    "green": "#ecfdf5",
    "gold":  "#fffbeb",
}
TEXT_DARK = "#1e293b"

LW, RW, H = 3.6, 3.9, 1.02          # box widths / height
LX, RX = 0.25, 4.55                 # column x
FS = 6.9                            # in-box font size

pipeline = [
    (9.05, "PubChemQC database\n3.2M molecules", "grey"),
    (7.75, "heteroatom + MW filter\n17,458 molecules", "grey"),
    (6.45, "B3LYP/6-31G*//PM6\nHOMO, LUMO, gap, dipole", "blue"),
    (5.15, "Scharber model\npredicted PCE", "blue"),
    (3.85, "RDKit SAScore\nPCE$_{\\mathrm{SAScore}}$ = PCE $-$ SAScore", "blue"),
    (2.55, "rank by PCE$_{\\mathrm{SAScore}}$\n7 molecules with score $>$ 0", "grey"),
]
gate = [
    (6.45, "raw viable set: 7\n(O$_2$, MgCO$_3$, cocrystal, azides)", "red"),
    (5.05, "chemical-validity filter\nsingle covalent $\\cdot$ aromatic ring\n$\\geq$ 6 conjugated atoms", "green"),
    (3.75, "chemically valid viable: 4\n17851, 20778, 4550, 1712", "green"),
    (2.45, "stability screen\nremove azide / nitroso / diazo", "green"),
    (1.15, "stable viable: 1\nmolecule 1712", "gold"),
]

fig, ax = plt.subplots(figsize=(fs.DOUBLE*0.63, 12.9/2.54))
ax.set_xlim(0, 8.6); ax.set_ylim(0.35, 10.75); ax.axis("off")
fig.patch.set_facecolor("white")

# ---- soft column backdrops (subtle grouping cue) ------------------------
backdrop_pad = 0.22
ax.add_patch(FancyBboxPatch(
    (LX - backdrop_pad, 2.35), LW + 2 * backdrop_pad, 7.9,
    boxstyle="round,pad=0,rounding_size=0.16",
    linewidth=0, facecolor="#f8fafc", zorder=0))
ax.add_patch(FancyBboxPatch(
    (RX - backdrop_pad, 0.95), RW + 2 * backdrop_pad, 6.7,
    boxstyle="round,pad=0,rounding_size=0.16",
    linewidth=0, facecolor="#f8fafc", zorder=0))


def shadow(x, y, w):
    """A faint offset duplicate of the box footprint, for gentle depth."""
    ax.add_patch(FancyBboxPatch(
        (x + 0.035, y - 0.045), w, H,
        boxstyle="round,pad=0.02,rounding_size=0.09",
        linewidth=0, facecolor="#0f172a", alpha=0.06, zorder=1))


def box(x, y, w, text, key, emphasize=False):
    ec, fc = INK[key], FILL[key]
    shadow(x, y, w)
    ax.add_patch(FancyBboxPatch(
        (x, y), w, H, boxstyle="round,pad=0.02,rounding_size=0.09",
        linewidth=1.6 if emphasize else 1.1,
        edgecolor=ec, facecolor=fc, zorder=2))
    ax.text(x + w / 2, y + H / 2, text, ha="center", va="center",
            fontsize=FS + (0.3 if emphasize else 0), color=TEXT_DARK,
            fontweight="bold" if emphasize else "normal", zorder=3)
    return x + w / 2, y, y + H


def arrow(x1, y1, x2, y2, color, cs=None):
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle="-|>,head_length=1.5,head_width=1",
        mutation_scale=6.5, linewidth=1.1, color=color,
        connectionstyle=cs, zorder=1.5, shrinkA=0, shrinkB=0))


# ---- pipeline column -----------------------------------------------------
pc = [box(LX, y, LW, t, k) for y, t, k in pipeline]
for a, b in zip(pc, pc[1:]):
    arrow(a[0], a[1] - 0.03, b[0], b[2] + 0.03, INK["grey"])

# ---- validity-gate column -------------------------------------------------
gc = [box(RX, y, RW, t, k, emphasize=(k == "gold")) for y, t, k in gate]
for a, b in zip(gc, gc[1:]):
    arrow(a[0], a[1] - 0.03, b[0], b[2] + 0.03, INK["green"])

# ranked output -> raw viable set (curved to skirt the columns' gap)
arrow(LX + LW + 0.03, 2.55 + H / 2, RX - 0.03, 6.45 + H / 2, INK["red"],
      cs="arc3,rad=-0.05")

# ---- section headers with a thin accent rule -----------------------------
def header(cx, y, text, color):
    ax.text(cx, y, text, fontsize=9, weight="bold", ha="center", color=color)
    ax.plot([cx - 1.15, cx + 1.15], [y - 0.22, y - 0.22],
            color=color, linewidth=1.6, solid_capstyle="round")

header(LX + LW / 2, 10.5, "Screening pipeline", INK["grey"])
header(RX + RW / 2, 7.9, "Validity gate (this work)", INK["green"])

fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"figures/figure1_workflow.{ext}", dpi=600, bbox_inches="tight")
print("wrote figures/figure1_workflow.pdf / .png")
