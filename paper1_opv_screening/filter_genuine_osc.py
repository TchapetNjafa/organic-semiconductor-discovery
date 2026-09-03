#!/usr/bin/env python3
"""
Chemical-validity filter for the PubChemQC OPV screen (rules R1-R6).

The raw PCE_SAScore ranking is contaminated by chemically implausible "molecules"
that pass the numeric gates but are not organic semiconductors (O2, MgCO3,
quinhydrone cocrystals, quaternary ammonium chlorides). This script applies a
transparent, reproducible set of chemistry rules to separate genuine conjugated
organic candidates from those artifacts, then re-ranks the survivors.

Rules (a molecule is REJECTED if any is true):
  R1  multi-fragment SMILES ('.')          -> salts / cocrystals / mixtures
  R2  contains a metal / metalloid         -> inorganic or organometallic salt
  R3  no ring at all                       -> OSCs need a (hetero)aromatic core
  R4  no aromatic ring                     -> need cyclic pi-conjugation
  R5  < 6 heavy atoms                      -> too small (O2, etc.)
  R6  < 6 sp2/aromatic (conjugated) atoms  -> insufficient pi-conjugation

Nothing here fabricates data: it only removes rows and re-ranks the rest.

INPUT.  outputs/recompute/screen_full.csv, written by scharber_recompute.py.
        Run that first. This script does NOT read the legacy property table:
        the efficiencies in that file came from a superseded implementation and
        must not re-enter the analysis.

The rules are also applied inside scharber_recompute.py, so screen_full.csv
already carries a validity_reject column. This script re-derives the rules from
the SMILES independently and asserts agreement, which makes it a check on the
pipeline rather than a duplicate of it.

Run:  python filter_genuine_osc.py
"""
import argparse
import pathlib
import sys

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = pathlib.Path(__file__).resolve().parent
DEFAULT_IN = HERE / "outputs" / "recompute" / "screen_full.csv"

# Metals / metalloids that mark a row as inorganic or organometallic.
METALS = {
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "Rb", "Sr", "Y", "Zr", "Nb",
    "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Cs", "Ba",
    "La", "Ce", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "B", "Si", "As", "Te",
}

# Which PCE_SAScore column to rank on. The four columns differ by absorber
# convention (strict / complementary) and fill factor (0.65 / Green); see the
# main text. The default is the complementary convention at FF = 0.65.
SCORE_COLUMNS = [
    "PCE_SAScore_compl_FF065",
    "PCE_SAScore_compl_FFgreen",
    "PCE_SAScore_strict_FF065",
    "PCE_SAScore_strict_FFgreen",
]


def reject_reason(smiles):
    """Return '' if the molecule is a plausible OSC, else the failing rule."""
    if not isinstance(smiles, str) or not smiles:
        return "unparsable"
    if "." in smiles:
        return "R1_multifragment"
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "unparsable"
    if {a.GetSymbol() for a in mol.GetAtoms()} & METALS:
        return "R2_metal"
    if mol.GetRingInfo().NumRings() == 0:
        return "R3_no_ring"
    if not any(a.GetIsAromatic() for a in mol.GetAtoms()):
        return "R4_no_aromatic_ring"
    if mol.GetNumHeavyAtoms() < 6:
        return "R5_too_small"
    conj = sum(
        1
        for a in mol.GetAtoms()
        if a.GetIsAromatic() or a.GetHybridization() == Chem.HybridizationType.SP2
    )
    if conj < 6:
        return "R6_low_conjugation"
    return ""


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--input", type=pathlib.Path, default=DEFAULT_IN,
                    help="screen_full.csv from scharber_recompute.py")
    ap.add_argument("--score", default=SCORE_COLUMNS[0], choices=SCORE_COLUMNS,
                    help="which PCE_SAScore convention to rank on")
    ap.add_argument("--out", type=pathlib.Path, default=HERE / "genuine_osc_ranked.csv")
    a = ap.parse_args()

    if not a.input.exists():
        sys.exit(f"missing {a.input}\nRun scharber_recompute.py first.")

    df = pd.read_csv(a.input, low_memory=False)
    df["reject"] = df["SMILES"].map(reject_reason)

    # Cross-check against the pipeline's own column: the two implementations must
    # agree on every one of the 17,458 rows.
    pipe = df["validity_reject"].fillna("").astype(str)
    disagree = df.loc[pipe != df["reject"], ["mol_id", "validity_reject", "reject"]]
    if len(disagree):
        print("PIPELINE DISAGREEMENT on", len(disagree), "rows:")
        print(disagree.head(10).to_string(index=False))
        sys.exit("validity rules do not match scharber_recompute.py -- fix before use")
    print(f"rule cross-check: {len(df)} rows agree with scharber_recompute.py")

    kept = df[df["reject"] == ""].sort_values(a.score, ascending=False)
    rejected = df[df["reject"] != ""]

    print(f"\nTotal molecules       : {len(df)}")
    print(f"Chemically valid      : {len(kept)}")
    print(f"Rejected as artifacts : {len(rejected)}")
    print("\nRejection breakdown:")
    print(rejected["reject"].value_counts().to_string())

    admitted = kept[kept[a.score].notna()]
    print(f"\n=== Valid molecules the physical gates admit, ranked by {a.score} ===")
    cols = ["mol_id", "formula", a.score, "SAScore", "gap_eV", "reactive_groups", "SMILES"]
    for _, r in admitted.head(15).iterrows():
        rg = r["reactive_groups"] if isinstance(r["reactive_groups"], str) and r["reactive_groups"] else "-"
        print(f"  {int(r['mol_id']):>6}  PS={r[a.score]:7.2f}  SAS={r['SAScore']:.2f}  "
              f"Eg={r['gap_eV']:.3f}  [{rg:<8}] {r['formula']:<14} {r['SMILES']}")

    print("\n=== Fate of the species that dominated the unfiltered ranking ===")
    for mid in [17851, 20778, 4550, 977, 11029, 1712, 7801, 24964, 9168]:
        row = df[df["mol_id"] == mid]
        if row.empty:
            print(f"  {mid:>6}  not in this dataset")
            continue
        row = row.iloc[0]
        verdict = "valid" if row["reject"] == "" else f"REJECTED ({row['reject']})"
        ps = row[a.score]
        ps_s = "gated out" if pd.isna(ps) else f"PS={ps:7.2f}"
        print(f"  {mid:>6} {str(row['formula']):<14} {ps_s:>14}  -> {verdict}")

    kept[["mol_id", "formula", a.score, "SAScore", "gap_eV", "HOMO_eV", "LUMO_eV",
          "reactive_groups", "stable", "SMILES"]].to_csv(a.out, index=False)
    print(f"\nWrote {len(kept)} chemically valid molecules -> {a.out.name}")


if __name__ == "__main__":
    sys.exit(main())
