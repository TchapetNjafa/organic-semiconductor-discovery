#!/usr/bin/env python3
"""
regen_fig6.py — heteroatom composition against frontier-orbital energy.

Replaces generate_fig6_heteroatom.py, which read a property table
(PCE_paper_GDB9.csv) that is no longer part of the deposit, and which marked
molecules 17851, 20778, 4550 and 1712 as "the four valid candidates" — a
labelling the corrected screen contradicts: 17851 and 20778 carry azides, 4550
carries a nitroso group, and 1712 is not admitted at all.

Everything here comes from outputs/recompute/screen_full.csv, so the r values in
the panels are the same numbers as tab_correlations.tex.

Emits figures/figure6_heteroatom_scatter.{pdf,png}.
Run from the Paper-1 folder:  python regen_fig6.py
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from rdkit import Chem
from rdkit import RDLogger

import figstyle as fs

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
REC = HERE / "outputs" / "recompute"
FIGDIR = HERE / "figures"
FIGDIR.mkdir(exist_ok=True)

df = pd.read_csv(REC / "screen_full.csv")
robust = json.loads((REC / "robustness.json").read_text())

INK = {"grey": "#a0aec0", "fit": "#c53030", "green": "#2f855a",
       "purple": "#6b46c1"}


def frac(smi, el):
    m = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
    if m is None:
        return np.nan
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    return 100 * syms.count(el) / len(syms) if syms else np.nan


df["fN"] = [frac(s, "N") for s in df.SMILES]
df["fO"] = [frac(s, "O") for s in df.SMILES]
d = df.dropna(subset=["fN", "fO", "HOMO_eV", "LUMO_eV"])

LEAD = robust["intersection_all_variants"][0] if robust["intersection_all_variants"] else None
OTHERS = [m for m in robust["union_all_variants"] if m != LEAD]

fig, axes = plt.subplots(1, 2, figsize=(fs.DOUBLE, fs.DOUBLE * 0.44))
panels = [("fN", "HOMO_eV", "nitrogen fraction (%)", "HOMO (eV)"),
          ("fO", "LUMO_eV", "oxygen fraction (%)", "LUMO (eV)")]

for ax, (xc, yc, xl, yl) in zip(axes, panels):
    ax.scatter(d[xc], d[yc], s=3, c=INK["grey"], alpha=0.35, linewidths=0,
               rasterized=True)
    r = float(np.corrcoef(d[xc], d[yc])[0, 1])
    b = np.polyfit(d[xc], d[yc], 1)
    xs = np.linspace(d[xc].min(), d[xc].max(), 50)
    ax.plot(xs, np.polyval(b, xs), color=INK["fit"], linewidth=1.6)
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    ax.text(0.04, 0.05, f"$r = {r:+.2f}$   $n = {len(d)}$",
            transform=ax.transAxes, fontsize=8, va="bottom", ha="left",
            bbox=dict(boxstyle="round", fc="white", ec="#cbd5e0"))
    print(f"  {xc} vs {yc}: r = {r:+.3f}  slope = {b[0]:+.4f} eV per % ")

# Mark the molecules the corrected screen actually admits, on both panels.
for ax, (xc, yc, *_ ) in zip(axes, panels):
    for mid in OTHERS:
        rr = df[df.mol_id == mid]
        if len(rr):
            ax.scatter(rr.iloc[0][xc], rr.iloc[0][yc], marker="s", s=26,
                       facecolor="none", edgecolor=INK["purple"], lw=1.2, zorder=5)
    if LEAD is not None:
        rr = df[df.mol_id == LEAD]
        if len(rr):
            ax.scatter(rr.iloc[0][xc], rr.iloc[0][yc], marker="o", s=44,
                       facecolor=INK["green"], edgecolor="black", lw=0.6, zorder=6)

axes[0].legend(handles=[
    Line2D([0], [0], marker="o", ls="", mfc=INK["green"], mec="black",
           label=f"robust lead ({LEAD})"),
    Line2D([0], [0], marker="s", ls="", mfc="none", mec=INK["purple"],
           label=f"other admitted, valid and stable ({len(OTHERS)})"),
    Line2D([0], [0], color=INK["fit"], lw=1.6, label="least-squares fit"),
], fontsize=6.4, loc="upper right", framealpha=0)

axes[0].text(0.03, 0.93, "(a)", transform=axes[0].transAxes, fontsize=10, weight="bold")
axes[1].text(0.03, 0.93, "(b)", transform=axes[1].transAxes, fontsize=10, weight="bold")
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(FIGDIR / f"figure6_heteroatom_scatter.{ext}", dpi=600,
                bbox_inches="tight")
plt.close(fig)
print(f"wrote figures/figure6_heteroatom_scatter.pdf / .png  (n = {len(d)})")
