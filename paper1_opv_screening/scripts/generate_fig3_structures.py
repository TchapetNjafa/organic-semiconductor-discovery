#!/usr/bin/env python3
"""
Figure 3: 2D structures of the four chemically valid viable candidates,
with reactive-group flags. Honest replacement for the old top-7 version.
Emits figures/figure3_structures.{pdf,png}. Data: viable_view_a_metric.csv.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import figstyle as fs
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image
import io


df = pd.read_csv("viable_view_a_metric.csv").sort_values("PCE_SAScore", ascending=False)

labels = {
    17851: "azide — unstable",
    20778: "azide — unstable",
    4550: "nitroso",
    1712: "stable",
}


def render(smiles, size=(430, 330)):
    mol = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(mol)
    d = rdMolDraw2D.MolDraw2DCairo(*size)
    opts = d.drawOptions()
    opts.padding = 0.12
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    return Image.open(io.BytesIO(d.GetDrawingText()))

fig, axes = plt.subplots(2, 2, figsize=(fs.DOUBLE, fs.DOUBLE*0.76))
for ax, (_, r) in zip(axes.ravel(), df.iterrows()):
    mid = int(r["mol_id"])
    ax.imshow(render(r["SMILES"]))
    ax.axis("off")
    flag = labels.get(mid, "")
    color = "#2f855a" if flag == "stable" else "#c53030"
    ax.set_title(
        f"Molecule {mid}  ({r['formula']})\n"
        f"PCE$_{{\\mathrm{{SAScore}}}}$ = {r['PCE_SAScore']:+.1f}   "
        f"SAScore = {r['sas1(%)']:.1f}",
        fontsize=8.5)
    ax.text(0.5, -0.06, flag, transform=ax.transAxes, ha="center", va="top",
            fontsize=8.5, color=color, weight="bold")

fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"figures/figure3_structures.{ext}", dpi=600, bbox_inches="tight")
print("wrote figures/figure3_structures.pdf / .png")
