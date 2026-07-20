#!/usr/bin/env python3
"""
Remaining Epic-1 figures for Paper 2:
  Fig 1 -- graphical abstract: two-axis schematic (docking / optical)
  Fig 2 -- workflow: Paper 1 screen -> 7 candidates -> docking + TD-DFT branches
  Fig 7 -- 2D structures of the 7 candidates, dipole MAGNITUDE labeled
           (no vector direction -- the source dataset only stores dipole
           magnitude, not x/y/z components, so a direction would be
           fabricated; magnitude-only annotation is what the data supports)
Style follows Paper 1's generate_fig1_workflow.py flowchart pattern
(matplotlib FancyBboxPatch/FancyArrowPatch, dependency-free) for series
consistency, plus RDKit for the 2D structure depictions in Fig 7.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from rdkit import Chem
from rdkit.Chem import Draw

import figstyle_p2 as fs

BASE = Path(__file__).resolve().parent

INK = {"grey": "#475569", "blue": "#1d4ed8", "green": "#047857", "gold": "#b45309", "red": "#b91c1c"}
FILL = {"grey": "#f1f5f9", "blue": "#eff6ff", "green": "#ecfdf5", "gold": "#fffbeb", "red": "#fef2f2"}
TEXT_DARK = "#1e293b"

SMILES = {
    "1712": "O=C1N=c2c(=C1CNc1ccc(cc1)S(=O)(=O)Nc1ccccn1)c1scnc1cc2",
    "9168": "c1ccc2c(c1)c1ccc3c(c1cc2)ccc1c3ccc2c1cccc2",
    "17574": "CC(=O)N(c1ccc(cc1)Nc1cc(c(c2c1C(=O)c1ccccc1C2=O)N)S(=O)(=O)O)C",
    "18506": "CCCCn1c(=O)c2c(c1=O)c(N)c1c(c2N)c(=O)c2c(c1=O)cccc2",
    "4550": "O=NC1=c2cc(N(O)O)c3c(c2=NC1=O)CCCC3",
    "17851": "CC[N-]c1nc(N[N+]#N)nc(n1)NC(C)C",
    "20778": "N#[N+]Nc1nc([N-]C(C)C)nc(n1)SC",
}
# Stable donors first, then reactive structures.
MOL_ORDER = ["1712", "9168", "17574", "18506", "4550", "17851", "20778"]
DIPOLE = {
    "1712": 2.56, "9168": 3.18, "17574": 5.18, "18506": 1.60,
    "4550": 9.74, "17851": 9.26, "20778": 7.42,
}


def _box(ax, x, y, w, h, text, color, fs_pt=7.0):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
        linewidth=1.1, edgecolor=INK[color], facecolor=FILL[color], zorder=2,
    ))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fs_pt, color=TEXT_DARK, zorder=3)


def _arrow(ax, x0, y0, x1, y1, color="grey"):
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=9,
        linewidth=1.1, color=INK[color], zorder=1,
    ))


def fig1_graphical_abstract():
    fig, ax = plt.subplots(figsize=(fs.SINGLE, fs.SINGLE * 0.85))
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

    _box(ax, 3.0, 8.0, 4.0, 1.3, "7 candidates\n(4 stable + 3 reactive)", "grey", 7.5)
    _arrow(ax, 4.3, 8.0, 2.3, 6.3, "green")
    _arrow(ax, 5.7, 8.0, 7.7, 6.3, "red")

    _box(ax, 0.3, 4.7, 4.0, 1.5, "4 stable donors\nbind in inhibitor band,\noptically insensitive", "green", 7.0)
    _box(ax, 5.7, 4.7, 4.0, 1.5, "3 reactive structures\nbind weaker; large shift\nin dark transitions", "red", 7.0)

    _arrow(ax, 2.3, 4.7, 2.3, 2.3, "green")
    _arrow(ax, 7.7, 4.7, 7.7, 2.3, "red")

    _box(ax, 0.3, 0.9, 4.0, 1.3, "Docking axis\n(molecular recognition)", "blue", 7.5)
    _box(ax, 5.7, 0.9, 4.0, 1.3, "Optical axis\n(solvatochromic shift)", "gold", 7.5)

    fs.save(fig, "fig1_graphical_abstract")


def fig2_workflow():
    fig, ax = plt.subplots(figsize=(fs.DOUBLE * 0.6, 12.5 / 2.54))
    ax.set_xlim(0, 8.6); ax.set_ylim(0.3, 10.6); ax.axis("off")

    _box(ax, 0.4, 9.0, 7.8, 1.1, "Paper 1: 17,458-molecule screen -> chemical-validity +\nstability gate; 4 admitted + 3 stable donors added = 7", "grey", 6.8)
    _arrow(ax, 2.3, 9.0, 2.3, 7.6, "blue")
    _arrow(ax, 6.3, 9.0, 6.3, 7.6, "gold")

    _box(ax, 0.4, 6.3, 4.0, 1.2, "Blind molecular docking\n(smina, seeded, 4 disease\ntargets)", "blue", 6.6)
    _box(ax, 4.6, 6.3, 3.6, 1.2, "TD-DFT solvatochromic shift\n(CPCM toluene/water,\nB3LYP + CAM-B3LYP)", "gold", 6.6)

    _arrow(ax, 2.3, 6.3, 2.3, 4.9, "blue")
    _arrow(ax, 6.3, 6.3, 6.3, 4.9, "gold")

    _box(ax, 0.4, 3.7, 4.0, 1.1, "Fig 3: affinity per target\n4 stable donors in inhibitor band", "blue", 6.6)
    _box(ax, 4.6, 3.7, 3.6, 1.1, "Fig 5: S1 shift per solvent\nlarge shift only in reactive set", "gold", 6.6)

    _arrow(ax, 2.3, 3.7, 4.3, 2.4, "grey")
    _arrow(ax, 6.3, 3.7, 4.3, 2.4, "grey")

    _box(ax, 1.3, 1.1, 6.0, 1.2, "Fig 6: two independent readouts,\nreported side by side (not merged)", "green", 7.0)

    fs.save(fig, "fig2_workflow")


def fig7_structures():
    # 7 candidates in a 2x4 grid (stable donors first, reactive last);
    # the 8th cell is blank.
    fig, axes = plt.subplots(2, 4, figsize=(fs.DOUBLE, fs.DOUBLE * 0.5))
    flat = axes.flatten()
    for ax, mol_id in zip(flat, MOL_ORDER):
        mol = Chem.MolFromSmiles(SMILES[mol_id])
        img = Draw.MolToImage(mol, size=(400, 400))
        ax.imshow(img)
        ax.axis("off")
        cls = "stable" if mol_id in ("1712", "9168", "17574", "18506") else "reactive"
        ax.set_title(f"{mol_id} ({cls})\n$\\mu$ = {DIPOLE[mol_id]:.2f} D", fontsize=8)
    for ax in flat[len(MOL_ORDER):]:
        ax.axis("off")
    fs.save(fig, "fig7_structures_dipole")


if __name__ == "__main__":
    (BASE / "figures").mkdir(exist_ok=True)
    fig1_graphical_abstract()
    fig2_workflow()
    fig7_structures()
