#!/usr/bin/env python3
"""
Generates Fig 3 (docking affinities), Fig 5 (solvatochromic shift, both
functionals), Fig 6 (dipole/docking/shift synthesis) for Paper 2, from the real
data files produced this session:
  docking_rerun/docking_results_clean.csv
  solvatochromic_tddft/solvatochromic_results_{B3LYP,CAM-B3LYP}.csv
Dipole moments are the seven verified values from Paper 1's dataset (hardcoded
here as a small, fixed, already-cited set -- not re-derived).
"""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import figstyle_p2 as fs

BASE = Path(__file__).resolve().parent
# Ordered stable donors first (green/blue), then reactive structures (red/gold/orange).
MOLECULES = ["1712", "9168", "17574", "18506", "4550", "17851", "20778"]
# Debye, Paper 1 dataset (dataset_pubchemqc_opv_17458.csv, col dipole_moment).
DIPOLE = {
    "1712": 2.56, "9168": 3.18, "17574": 5.18, "18506": 1.60,
    "4550": 9.74, "17851": 9.26, "20778": 7.42,
}
TARGETS = ["HIV1protease", "Hsp90", "Neurodegenerative_1SYH", "COVID19_6Y2F"]
TARGET_LABELS = {
    "HIV1protease": "HIV-1 protease",
    "Hsp90": "Hsp90",
    "Neurodegenerative_1SYH": "Neurodeg.\n(1SYH)",
    "COVID19_6Y2F": "SARS-CoV-2\nMpro (6Y2F)",
}


def load_docking():
    path = BASE / "docking_rerun" / "docking_results_clean.csv"
    data = {mol: {} for mol in MOLECULES}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            data[row["mol_id"]][row["target"]] = float(row["affinity_kcal_mol"])
    return data


def load_solvatochromic(functional):
    path = BASE / "solvatochromic_tddft" / f"solvatochromic_results_{functional}.csv"
    s1 = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            if int(row["state"]) == 1:
                s1.setdefault(row["mol_id"], {})[row["solvent"]] = float(row["energy_eV"])
    return s1


def fig3_docking():
    data = load_docking()
    x = np.arange(len(TARGETS))
    n = len(MOLECULES)
    width = 0.8 / n
    fig, ax = plt.subplots(figsize=(fs.DOUBLE, fs.DOUBLE * 0.42))
    for i, mol in enumerate(MOLECULES):
        vals = [data[mol][t] for t in TARGETS]
        offset = (i - (n - 1) / 2) * width
        ax.bar(x + offset, vals, width, label=mol, color=fs.MOL_COLORS[mol])
    ax.set_xticks(x)
    ax.set_xticklabels([TARGET_LABELS[t] for t in TARGETS])
    ax.set_ylabel("Docking affinity (kcal/mol)")
    ax.invert_yaxis()
    ax.legend(title="Candidate", ncol=n, loc="upper center", bbox_to_anchor=(0.5, 1.18))
    fs.save(fig, "fig3_docking_affinities")


def fig5_solvatochromic():
    fig, axes = plt.subplots(1, 2, figsize=(fs.DOUBLE, fs.DOUBLE * 0.42), sharey=False)
    for ax, functional, title in zip(axes, ["B3LYP", "CAM-B3LYP"], ["B3LYP/6-31G*", "CAM-B3LYP/6-31G*"]):
        s1 = load_solvatochromic(functional)
        x = np.arange(len(MOLECULES))
        tol = [s1[m]["Toluene"] for m in MOLECULES]
        wat = [s1[m]["Water"] for m in MOLECULES]
        width = 0.35
        ax.bar(x - width / 2, tol, width, label="Toluene", color=fs.COLORS["light"], edgecolor=fs.COLORS["grey"])
        ax.bar(x + width / 2, wat, width, label="Water", color=fs.COLORS["primary"])
        ax.set_xticks(x)
        ax.set_xticklabels(MOLECULES)
        ax.set_ylabel("S1 energy (eV)")
        ax.set_title(title)
    axes[0].legend(loc="upper left")
    fs.save(fig, "fig5_solvatochromic_shift")


def fig6_synthesis():
    docking = load_docking()
    best_dock = {m: min(docking[m].values()) for m in MOLECULES}  # most negative = strongest
    shift_b3lyp = load_solvatochromic("B3LYP")
    shift_camb3lyp = load_solvatochromic("CAM-B3LYP")

    def shift_mag(s1, mol):
        return abs(s1[mol]["Water"] - s1[mol]["Toluene"])

    fig, ax1 = plt.subplots(figsize=(fs.SINGLE, fs.SINGLE * 0.9))
    dip = [DIPOLE[m] for m in MOLECULES]
    dock = [best_dock[m] for m in MOLECULES]
    ax1.scatter(dip, dock, color=fs.COLORS["primary"], s=40, label="Docking (best target)", zorder=3)
    for m, x, y in zip(MOLECULES, dip, dock):
        ax1.annotate(m, (x, y), textcoords="offset points", xytext=(4, 4), fontsize=7)
    ax1.set_xlabel("Dipole moment (D)")
    ax1.set_ylabel("Best docking affinity (kcal/mol)", color=fs.COLORS["primary"])
    ax1.invert_yaxis()

    ax2 = ax1.twinx()
    shift_cb = [shift_mag(shift_camb3lyp, m) for m in MOLECULES]
    shift_b3 = [shift_mag(shift_b3lyp, m) for m in MOLECULES]
    ax2.scatter(dip, shift_cb, color=fs.COLORS["artifact"], marker="^", s=40, label="Shift (CAM-B3LYP)", zorder=3)
    ax2.scatter(dip, shift_b3, color=fs.COLORS["artifact"], marker="s", s=25, facecolors="none", label="Shift (B3LYP)", zorder=2)
    ax2.set_ylabel("|Solvatochromic shift| (eV)", color=fs.COLORS["artifact"])

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper center", fontsize=6.5, bbox_to_anchor=(0.5, 1.22), ncol=1)
    fs.save(fig, "fig6_dipole_switch_synthesis")


if __name__ == "__main__":
    (BASE / "figures").mkdir(exist_ok=True)
    fig3_docking()
    fig5_solvatochromic()
    fig6_synthesis()
