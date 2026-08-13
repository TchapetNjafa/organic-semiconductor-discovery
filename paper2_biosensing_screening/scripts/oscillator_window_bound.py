#!/usr/bin/env python3
"""Oscillator-strength ceiling of the eight-root TD-DFT window.

Reproduces the bound reported in the Discussion (subsection "Environment
sensitivity and optical accessibility sit on different molecules") and recorded
in analysis-ledger entries L-2026-08-12-15 and L-2026-08-12-16.

The manuscript claim is deliberately narrow, and this script is written to show
exactly how narrow. At CAM-B3LYP/6-31G* no state of 17851 or 20778 above
400 nm reaches f = 0.04 in either solvent. That bound does NOT hold for 4550,
and it does NOT hold at other functional/basis combinations:

    CAM-B3LYP/6-31G*   4550 water    f = 0.617 at 392 nm  (below 400 nm)
    CAM-B3LYP/6-31+G*  4550 water    f = 0.164 at 415 nm  (BREACH)
    B3LYP/6-31G*       4550 toluene  f = 0.497 at 429 nm  (BREACH)
    B3LYP/6-31+G*      4550 water    f = 0.409 at 433 nm  (BREACH)

The script prints every breach so the scope of the manuscript sentence can be
checked rather than trusted. An earlier draft of that sentence claimed the bound
for all three molecules in "both solvents" without naming a method, which is
false; see L-2026-08-12-16.

The window ceiling itself is reported per molecule and solvent, because the
eight-root cutoff lands at a different energy for each and no single uniform
bound exists across the set.

Usage:
    python3 oscillator_window_bound.py [--data-dir DIR]
                                       [--wavelength-cut NM]
                                       [--f-threshold F]

Requires: pandas.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

# Molecules carrying the large solvatochromic shifts (the "reactive trio").
LARGE_SHIFT = (4550, 17851, 20778)

# Molecules the bound is claimed for in the manuscript.
BOUNDED = (17851, 20778)

DEFAULT_WAVELENGTH_CUT_NM = 400.0
DEFAULT_F_THRESHOLD = 0.04
EXPECTED_STATES_PER_GROUP = 8

FILE_LABELS = {
    "solvatochromic_results_B3LYP.csv": "B3LYP/6-31G*",
    "solvatochromic_results_B3LYP_diffuse.csv": "B3LYP/6-31+G*",
    "solvatochromic_results_CAM-B3LYP.csv": "CAM-B3LYP/6-31G*",
    "solvatochromic_results_CAM-B3LYP_diffuse.csv": "CAM-B3LYP/6-31+G*",
    "solvatochromic_results_gap_candidates_B3LYP.csv": "B3LYP/6-31G* (gap)",
    "solvatochromic_results_gap_candidates_CAM-B3LYP.csv": "CAM-B3LYP/6-31G* (gap)",
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="directory holding the solvatochromic_results_*.csv files",
    )
    parser.add_argument(
        "--wavelength-cut",
        type=float,
        default=DEFAULT_WAVELENGTH_CUT_NM,
        help="report states with wavelength longer than this, in nm",
    )
    parser.add_argument(
        "--f-threshold",
        type=float,
        default=DEFAULT_F_THRESHOLD,
        help="oscillator strength above which a state counts as a breach",
    )
    args = parser.parse_args()

    breaches = []

    for filename, label in sorted(FILE_LABELS.items()):
        path = args.data_dir / filename
        if not path.exists():
            print(f"skipping absent file: {filename}", file=sys.stderr)
            continue
        frame = pd.read_csv(path)

        print(f"=== {label}  ({filename})")
        for keys, group in frame.groupby(["mol_id", "solvent"]):
            mol_id = int(keys[0])
            solvent = str(keys[1])
            if mol_id not in LARGE_SHIFT:
                continue
            if len(group) != EXPECTED_STATES_PER_GROUP:
                print(
                    f"    WARNING {mol_id}/{solvent}: {len(group)} states, "
                    f"expected {EXPECTED_STATES_PER_GROUP}"
                )

            ceiling = float(group.energy_eV.max())
            long_wave = group[group.wavelength_nm > args.wavelength_cut]
            if long_wave.empty:
                print(
                    f"    {mol_id:>6} {solvent:<8} ceiling "
                    f"{ceiling:5.2f} eV   no states above "
                    f"{args.wavelength_cut:.0f} nm"
                )
                continue

            peak = long_wave.loc[long_wave.oscillator_strength.idxmax()]
            peak_f = float(peak.oscillator_strength)
            peak_nm = float(peak.wavelength_nm)
            flag = ""
            if peak_f > args.f_threshold:
                flag = "  <-- BREACH"
                breaches.append(
                    f"{label} {mol_id}/{solvent}: "
                    f"f={peak_f:.4f} at {peak_nm:.0f} nm"
                )
            print(
                f"    {mol_id:>6} {solvent:<8} ceiling {ceiling:5.2f} eV   "
                f"max f above {args.wavelength_cut:.0f} nm = "
                f"{peak_f:.4f} at {peak_nm:5.0f} nm{flag}"
            )
        print()

    print(
        f"--- summary: {len(breaches)} state(s) above "
        f"{args.wavelength_cut:.0f} nm exceed f = {args.f_threshold}"
    )
    for line in breaches:
        print(f"    {line}")
    print(
        "\nThe manuscript claims the bound only for "
        f"{' and '.join(str(m) for m in BOUNDED)} at CAM-B3LYP/6-31G*. "
        "Every breach listed above lies outside that scope; a breach appearing "
        "inside it would falsify the manuscript sentence."
    )


if __name__ == "__main__":
    main()

# DATA PROVENANCE (deposited copy).
# This repository holds code and geometries; the full CSV datasets are archived
# on Zenodo under concept DOI 10.5281/zenodo.18201812. The diffuse-basis and
# gap-candidate result files this script reads are Zenodo-only. Download them
# and pass their location with --data-dir, e.g.
#     python3 oscillator_window_bound.py --data-dir /path/to/zenodo/data
