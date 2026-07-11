#!/usr/bin/env python3
"""
Figure 8: sensitivity of the viable-set size to the SAScore weight.
Honest note: these counts are PRE chemical-validity filter, so they still
include artifacts (e.g. MgCO3, id 11029). Data: ZENODO_UPLOAD/pce_sascore_sensitivity_results.csv.
Emits figures/figure8_sensitivity.{pdf,png}.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import figstyle as fs
import pandas as pd


df = pd.read_csv("../../data/pce_sascore_sensitivity_results.csv")
# ratio of SA weight to PCE weight controls how hard synthesis is penalised
df["w"] = df["SA_weight"] / df["PCE_weight"]
df = df.sort_values("w")

fig, ax = plt.subplots(figsize=(fs.SINGLE, fs.SINGLE*0.84))
ax.plot(df["w"], df["N_candidates"], "o-", color="#2b6cb0", linewidth=1.6, markersize=6)
for _, r in df.iterrows():
    ax.annotate(f"{int(r['N_candidates'])}", (r["w"], r["N_candidates"]),
                textcoords="offset points", xytext=(0, 7), ha="center", fontsize=8)
ax.set_xlabel("SAScore weight / PCE weight  ($w$)")
ax.set_ylabel("molecules with PCE$_{\\mathrm{SAScore}} > 0$")
ax.set_ylim(0, 5)
ax.axhline(1, color="#c53030", linestyle="--", linewidth=1,
           label="stable viable after validity + stability screen (1)")
ax.text(1.55, 2.55,
        "counts are pre-validity-filter\n(still include O$_2$, MgCO$_3$, cocrystal)",
        ha="left", va="center", fontsize=7.0, color="#4a5568")
ax.legend(fontsize=7.0, loc="lower left")
ax.grid(alpha=0.25)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"figures/figure8_sensitivity.{ext}", dpi=600, bbox_inches="tight")
print("wrote figures/figure8_sensitivity.pdf / .png")
