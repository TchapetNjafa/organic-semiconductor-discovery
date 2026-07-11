#!/usr/bin/env python3
"""
Figure 7: corrected four-point hole reorganization energy for the subset with
optimized charged-state geometries, across three functionals. The azido-triazines
17851 and 20778 give NEGATIVE lambda_hole (unphysical -> instability under
ionization). Data: NTO_and_related_study/analysis/reorganization_energies_corrected.csv.
Emits figures/figure7_reorganization_energies.{pdf,png}.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import figstyle as fs
import numpy as np
import pandas as pd


df = pd.read_csv("../../data/analysis/reorganization_energies_corrected.csv")
df = df[df["status"] == "success"]
funcs = ["b3lyp", "camb3lyp", "wb97x-d3"]
fcolor = {"b3lyp": "#2b6cb0", "camb3lyp": "#2f855a", "wb97x-d3": "#d69e2e"}
flabel = {"b3lyp": "B3LYP", "camb3lyp": "CAM-B3LYP", "wb97x-d3": "$\\omega$B97X-D3"}

# Show the two top-ranked azido-triazine candidates (negative, unphysical) against
# two well-behaved neutral organic controls (positive, physical). O2 (977) and the
# rejected chloride salt 19531 are excluded to avoid mixing artifacts into the panel.
mols = [17851, 20778, 14380, 18985]
mtag = {17851: "17851\n(azide)", 20778: "20778\n(azide)",
        14380: "14380\n(control)", 18985: "18985\n(control)"}

x = np.arange(len(mols))
w = 0.26
fig, ax = plt.subplots(figsize=(fs.SINGLE, fs.SINGLE*0.90))
for k, f in enumerate(funcs):
    vals = [df[(df.mol_id == m) & (df.functional == f)]["lambda_hole_eV"].values for m in mols]
    vals = [v[0] if len(v) else np.nan for v in vals]
    ax.bar(x + (k - 1) * w, vals, w, color=fcolor[f], label=flabel[f], edgecolor="black", linewidth=0.3)

ax.axhline(0, color="black", linewidth=0.8)
ax.axhspan(-1.2, 0, color="#c53030", alpha=0.07)
ax.text(0.5, -0.75, "negative $\\lambda_{\\mathrm{hole}}$ (unphysical)",
        ha="center", va="center", fontsize=7.2, color="#c53030")
ax.set_xticks(x)
ax.set_xticklabels([mtag[m] for m in mols], fontsize=7.5)
ax.set_ylabel("$\\lambda_{\\mathrm{hole}}$ (eV)")
ax.set_xlabel("molecule")
ax.legend(fontsize=7, ncol=3, loc="upper center", frameon=False)
ax.set_ylim(-1.2, 1.2)
ax.grid(axis="y", alpha=0.25)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"figures/figure7_reorganization_energies.{ext}", dpi=600, bbox_inches="tight")
print("wrote figures/figure7_reorganization_energies.pdf / .png")
