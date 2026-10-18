#!/usr/bin/env python3
"""Regenerate Figure 7: hole reorganization energy, all six molecules (T0-5).

WHY THIS REPLACES generate_fig7_reorg.py
----------------------------------------
The submitted figure showed 17851 and 20778 (azides, negative lambda_hole) against
two positive controls and explicitly EXCLUDED 19531, with the comment "to avoid
mixing artifacts into the panel". 19531's lambda_hole is -3.44 to -4.21 eV, an
order of magnitude beyond either azide.

The real diagnosis (verified, ledger L-2026-09-01-22): 19531 is
`CCCOc1cc(C)c(c(c1)C)C(=O)OCC([NH+]1CCCC1)C.[Cl-]` -- a pyrrolidinium ester
CHLORIDE ION PAIR, C19H30ClNO3, which the validity filter rejects as
R1_multifragment. Removing an electron from an ion pair transfers charge between
fragments rather than ionising a neutral molecule, so the four-point formula does
not describe a reorganization at all. Its term1 (-5.64 eV) is comparable to the
other molecules' but term2 is only +2.20 eV instead of +5.6 to +7.9 eV, and the
sum goes badly negative. The `anion_neutral_failed` flag on its wB97X-D3 row is a
SEPARATE issue and cannot explain a negative lambda_HOLE, which is a cation
quantity -- an error in my own earlier reading, corrected here.

Excluding the strongest instance of the effect, on the grounds that it does not
fit the explanation offered for the effect, is not a presentational choice. The
figure now shows all six molecules on a broken axis, states why 19531 is
pathological, and lets the reader see that negative lambda_hole is not
azide-specific.

Emits figures/figure7_reorganization_energies.{pdf,png}.
Run from the Paper-1 folder:  python regen_fig7.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import figstyle as fs

HERE = Path(__file__).resolve().parent
# In the working tree the corrected table lives in the sibling NTO study; in the
# published deposit a copy sits under data/. Both are accepted.
SRC = next((p for p in (HERE / "data" / "reorganization_energies_corrected.csv",
                        HERE.parent / "NTO_and_related_study" / "analysis" /
                        "reorganization_energies_corrected.csv") if p.exists()),
           HERE / "data" / "reorganization_energies_corrected.csv")
FIGDIR = HERE / "figures"

df = pd.read_csv(SRC)
FUNCS = ["b3lyp", "camb3lyp", "wb97x-d3"]
FCOLOR = {"b3lyp": "#2b6cb0", "camb3lyp": "#2f855a", "wb97x-d3": "#d69e2e"}
FLABEL = {"b3lyp": "B3LYP", "camb3lyp": "CAM-B3LYP", "wb97x-d3": "$\\omega$B97X-D3"}

piv = {f: df[df.functional == f].set_index("mol_id") for f in FUNCS}
ids = sorted(set(df.mol_id))
# order by B3LYP lambda_hole so the outlier sits at one end and reads as an outlier
order = sorted(ids, key=lambda m: piv["b3lyp"].loc[m, "lambda_hole_eV"])

TAG = {17851: "17851\nazide", 20778: "20778\nazide",
       19531: "19531\nion pair", 977: "977\n$\\mathrm{O_2}$",
       14380: "14380\ncontrol", 18985: "18985\ncontrol"}

x = np.arange(len(order))
w = 0.26

# broken axis: 19531 sits near -4 eV, everything else within +/-1 eV
fig, (top, bot) = plt.subplots(
    2, 1, sharex=True, figsize=(fs.SINGLE, fs.SINGLE * 1.02),
    gridspec_kw=dict(height_ratios=[2.6, 1.0], hspace=0.10))

for ax in (top, bot):
    for k, f in enumerate(FUNCS):
        vals = [piv[f].loc[m, "lambda_hole_eV"] for m in order]
        ax.bar(x + (k - 1) * w, vals, w, color=FCOLOR[f],
               label=FLABEL[f] if ax is top else None,
               edgecolor="black", linewidth=0.3)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.grid(axis="y", alpha=0.25)

top.set_ylim(-1.15, 1.25)
bot.set_ylim(-4.55, -3.15)
top.axhspan(-1.15, 0, color="#c53030", alpha=0.07)
bot.axhspan(-4.55, 0, color="#c53030", alpha=0.07)

# hide the shared spine and draw the break marks
top.spines["bottom"].set_visible(False)
bot.spines["top"].set_visible(False)
top.tick_params(bottom=False)
kw = dict(marker=[(-1, -0.6), (1, 0.6)], markersize=6, linestyle="none",
          color="black", mec="black", mew=0.8, clip_on=False)
top.plot([0, 1], [0, 0], transform=top.transAxes, **kw)
bot.plot([0, 1], [1, 1], transform=bot.transAxes, **kw)

# mark the real pathology: 19531 is an ion pair, not a neutral molecule
i19531 = order.index(19531)
bot.annotate("chloride ion pair;\nfour-point formula\nnot applicable",
             xy=(i19531 + w, piv["camb3lyp"].loc[19531, "lambda_hole_eV"]),
             xytext=(i19531 + 1.30, -3.62), fontsize=6.0, color="#c53030",
             ha="left", va="center",
             arrowprops=dict(arrowstyle="-", color="#c53030", lw=0.7))

top.text(len(order) - 1.4, -0.85, "negative $\\lambda_{\\mathrm{hole}}$\n(unphysical)",
         ha="center", va="center", fontsize=6.6, color="#c53030")
top.legend(fontsize=6.6, ncol=3, loc="upper left", frameon=False)
bot.set_xticks(x)
bot.set_xticklabels([TAG[m] for m in order], fontsize=6.8)
bot.set_xlabel("molecule")
fig.supylabel("$\\lambda_{\\mathrm{hole}}$ (eV)", fontsize=8.5, x=0.005)

fig.subplots_adjust(left=0.19, right=0.98, top=0.97, bottom=0.16)
for ext in ("pdf", "png"):
    fig.savefig(FIGDIR / f"figure7_reorganization_energies.{ext}", dpi=600)
plt.close(fig)

neg = [m for m in order if piv["b3lyp"].loc[m, "lambda_hole_eV"] < 0]
print(f"  figure7_reorg         {len(order)} molecules (was 4), negative lambda: {neg}")
print(f"                        19531 B3LYP {piv['b3lyp'].loc[19531, 'lambda_hole_eV']:+.3f} eV, "
      f"wB97X-D3 {piv['wb97x-d3'].loc[19531, 'lambda_hole_eV']:+.3f} eV "
      f"({piv['wb97x-d3'].loc[19531, 'status']})")
