#!/usr/bin/env python3
"""Spearman rank correlation between ground-state dipole moment and
solvatochromic shift magnitude, over the ten-candidate set.

Reproduces the values reported in the Discussion (section "Polarity orders
binding not at all, and optical response only in part") and recorded in
analysis-ledger entry L-2026-08-12-17:

    B3LYP      n=10  rho=0.758  permutation p=0.0076
    CAM-B3LYP  n=10  rho=0.636  permutation p=0.0277

Shift magnitude is |E(S1, water) - E(S1, toluene)|. Dipole moments are the
B3LYP/6-31G* ground-state values tabulated in the manuscript and its SI.

NOTE ON THE ASSERTION BELOW. An earlier ad-hoc version of this analysis silently
dropped two molecules (20549 and 23424) because its dipole lookup table had no
entry for them, and reported an n=8 correlation (rho=0.69/0.62) as though it
covered all ten candidates. The assertion on N_CANDIDATES exists so that a
missing molecule fails loudly instead of producing a plausible wrong number.

Usage:
    python3 correlation_dipole_shift.py [--data-dir DIR]

Requires: pandas, scipy.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import permutation_test, spearmanr

# Ground-state dipole moments, B3LYP/6-31G*, in debye.
# Source: manuscript Results (sections/results.tex) and SI section S1.
DIPOLE_DEBYE = {
    18506: 1.60,
    5081: 2.42,
    1712: 2.56,
    20549: 2.96,
    9168: 3.18,
    17574: 5.18,
    20778: 7.42,
    23424: 8.17,
    17851: 9.26,
    4550: 9.74,
}

N_CANDIDATES = 10
FUNCTIONALS = ("B3LYP", "CAM-B3LYP")
N_RESAMPLES = 20000
RANDOM_STATE = 0


def shift_table(data_dir: Path, functional: str) -> pd.DataFrame:
    """Return one row per molecule: mol_id, dipole (D), |shift| (eV)."""
    sources = [
        data_dir / f"solvatochromic_results_{functional}.csv",
        data_dir / f"solvatochromic_results_gap_candidates_{functional}.csv",
    ]
    rows = []
    for path in sources:
        if not path.exists():
            sys.exit(f"missing input: {path}")
        frame = pd.read_csv(path)
        for mol_id, group in frame.groupby("mol_id"):
            first = group[group.state == 1]
            toluene = first[first.solvent == "Toluene"]
            water = first[first.solvent == "Water"]
            if toluene.empty or water.empty:
                continue
            if mol_id not in DIPOLE_DEBYE:
                sys.exit(
                    f"molecule {mol_id} has no tabulated dipole moment; "
                    "add it to DIPOLE_DEBYE rather than letting it be dropped"
                )
            shift = abs(
                water.energy_eV.iloc[0] - toluene.energy_eV.iloc[0]
            )
            rows.append((int(mol_id), DIPOLE_DEBYE[mol_id], shift))

    table = pd.DataFrame(rows, columns=["mol_id", "dipole_D", "shift_eV"])
    table = table.sort_values("dipole_D").reset_index(drop=True)

    assert len(table) == N_CANDIDATES, (
        f"{functional}: expected {N_CANDIDATES} molecules, got {len(table)} "
        f"({sorted(table.mol_id)}). Refusing to report a correlation over an "
        "incomplete set."
    )
    return table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="directory holding the solvatochromic_results_*.csv files",
    )
    args = parser.parse_args()

    for functional in FUNCTIONALS:
        table = shift_table(args.data_dir, functional)
        rho, asymptotic_p = spearmanr(table.dipole_D, table.shift_eV)
        result = permutation_test(
            (table.dipole_D.values, table.shift_eV.values),
            lambda x, y: spearmanr(x, y).statistic,
            permutation_type="pairings",
            alternative="greater",
            n_resamples=N_RESAMPLES,
            random_state=RANDOM_STATE,
        )
        print(
            f"{functional}: n={len(table)} rho={rho:.3f} "
            f"asymptotic_p={asymptotic_p:.4f} permutation_p={result.pvalue:.4f}"
        )
        for row in table.itertuples():
            print(
                f"    {row.mol_id:>6}  {row.dipole_D:5.2f} D  "
                f"{row.shift_eV:6.3f} eV"
            )
        print()


if __name__ == "__main__":
    main()
