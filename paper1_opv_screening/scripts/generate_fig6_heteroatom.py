#!/usr/bin/env python3
"""
Figure 6: heteroatom composition vs frontier-orbital energy over the FULL
dataset (n=17,457) -- honest replacement for the Day-2 version, which computed
correlations on only 17 molecules and so reported spuriously strong trends.
(a) N fraction vs HOMO; (b) O fraction vs LUMO. Both correlations are weak.
Emits figures/figure6_heteroatom_scatter.{pdf,png}.
"""
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
import figstyle as fs
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")


df = pd.read_csv("../../data/PCE_paper_GDB9.csv")


def frac(smi, el):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return np.nan
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    return 100 * syms.count(el) / len(syms) if syms else np.nan


df["N%"] = df["SMILES"].map(lambda s: frac(s, "N"))
df["O%"] = df["SMILES"].map(lambda s: frac(s, "O"))
d = df.dropna(subset=["N%", "O%", "HOMO(eV)", "LUMO(eV)"])

fig, axes = plt.subplots(1, 2, figsize=(fs.DOUBLE, fs.DOUBLE*0.44))
panels = [("N%", "HOMO(eV)", "nitrogen fraction (%)", "HOMO (eV)"),
          ("O%", "LUMO(eV)", "oxygen fraction (%)", "LUMO (eV)")]
for ax, (xc, yc, xl, yl) in zip(axes, panels):
    ax.scatter(d[xc], d[yc], s=3, c="#a0aec0", alpha=0.35, linewidths=0, rasterized=True)
    r = d[xc].corr(d[yc])
    b = np.polyfit(d[xc], d[yc], 1)
    xs = np.linspace(d[xc].min(), d[xc].max(), 50)
    ax.plot(xs, np.polyval(b, xs), color="#c53030", linewidth=1.6)
    ax.set_xlabel(xl); ax.set_ylabel(yl)
    ax.text(0.04, 0.05, f"$r = {r:+.2f}$", transform=ax.transAxes, fontsize=9,
            va="bottom", ha="left", bbox=dict(boxstyle="round", fc="white", ec="#cbd5e0"))
# mark the four valid candidates on panel (a)
for mid in (17851, 20778, 4550, 1712):
    r0 = df[df.mol_id == mid].iloc[0]
    axes[0].scatter(r0["N%"], r0["HOMO(eV)"], marker="o", s=30, facecolor="#2f855a",
                    edgecolor="black", linewidths=0.5, zorder=5)
axes[0].text(0.03, 0.93, "(a)", transform=axes[0].transAxes, fontsize=10, weight="bold")
axes[1].text(0.03, 0.93, "(b)", transform=axes[1].transAxes, fontsize=10, weight="bold")
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"figures/figure6_heteroatom_scatter.{ext}", dpi=600, bbox_inches="tight")
print("wrote figures/figure6_heteroatom_scatter.pdf / .png")
