#!/usr/bin/env python3
"""
Figure 2: PCE vs SAScore landscape (honest version).
Every molecule plotted; the diagonal PCE = SAScore is the PCE_SAScore = 0 line.
The four chemically valid viable candidates are highlighted; the three artifacts
that also cross the threshold (O2, MgCO3, quinhydrone cocrystal) are marked
separately. Data: ZENODO_UPLOAD/PCE_paper_GDB9.csv.
Emits figures/figure2_pce_sascore_scatter.{pdf,png}.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import figstyle as fs
import pandas as pd


df = pd.read_csv("../../data/PCE_paper_GDB9.csv")
df["PCE"] = df[["pce_pcbm(%)", "pce_pcdtbt(%)"]].max(axis=1)
df["SAS"] = df["sas1(%)"]

valid = {17851: "17851", 20778: "20778", 4550: "4550", 1712: "1712"}
artifact = {977: "O$_2$", 11029: "MgCO$_3$", 7801: "cocrystal"}

fig, ax = plt.subplots(figsize=(fs.SINGLE, fs.SINGLE*0.92))
ax.scatter(df["SAS"], df["PCE"], s=4, c="#cbd5e0", alpha=0.5, linewidths=0, rasterized=True)

# viability diagonal PCE = SAScore
lim = [0, 10]
ax.plot(lim, lim, color="#2d3748", linestyle="--", linewidth=1,
        label="PCE$_{\\mathrm{SAScore}} = 0$")

# per-label offsets so the cocrystal marker (near molecule 1712) stays legible
art_off = {977: (6, -2), 11029: (6, -2), 7801: (-4, -13)}
art_ha = {977: "left", 11029: "left", 7801: "right"}
for mid, name in artifact.items():
    r = df[df["mol_id"] == mid].iloc[0]
    ax.scatter(r["SAS"], r["PCE"], marker="x", s=55, c="#c53030", linewidths=1.8, zorder=5)
    ax.annotate(name, (r["SAS"], r["PCE"]), textcoords="offset points",
                xytext=art_off[mid], ha=art_ha[mid], fontsize=7.5, color="#c53030")

for mid, name in valid.items():
    r = df[df["mol_id"] == mid].iloc[0]
    ax.scatter(r["SAS"], r["PCE"], marker="o", s=42, facecolor="#2f855a",
               edgecolor="black", linewidths=0.6, zorder=6)
    ax.annotate(name, (r["SAS"], r["PCE"]), textcoords="offset points",
                xytext=(6, 3), fontsize=7.5, color="#2f855a", weight="bold")

ax.set_xlabel("SAScore")
ax.set_ylabel("predicted PCE (%)")
ax.set_xlim(0.5, 10)
ax.set_ylim(-1, 40)
# legend proxies
from matplotlib.lines import Line2D
handles = [
    Line2D([0], [0], ls="--", color="#2d3748", label="PCE$_{\\mathrm{SAScore}} = 0$"),
    Line2D([0], [0], marker="o", ls="", mfc="#2f855a", mec="black", label="chemically valid viable (4)"),
    Line2D([0], [0], marker="x", ls="", color="#c53030", label="artifact above threshold (3)"),
]
ax.legend(handles=handles, fontsize=7.2, loc="upper left")
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"figures/figure2_pce_sascore_scatter.{ext}", dpi=600, bbox_inches="tight")
print("wrote figures/figure2_pce_sascore_scatter.pdf / .png")
