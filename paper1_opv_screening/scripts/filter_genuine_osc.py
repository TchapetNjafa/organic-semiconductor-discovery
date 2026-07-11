#!/usr/bin/env python3
"""
Chemical sanity filter for the PubChemQC OPV screen.

The raw PCE_SAScore ranking is contaminated by chemically implausible
"molecules" that pass the numeric filter but are not organic semiconductors
(e.g. O2, MgCO3 salts, quinhydrone cocrystals). This script applies a
transparent, reproducible set of chemistry rules to separate genuine
conjugated organic candidates from these artifacts, then re-ranks the
survivors by PCE_SAScore.

Rules (a molecule is REJECTED if any is true):
  R1  multi-fragment SMILES ('.')            -> salts / cocrystals / mixtures
  R2  contains a metal / metalloid           -> inorganic or organometallic salt
  R3  no ring at all                         -> OSCs require a (hetero)aromatic core
  R4  no aromatic ring                        -> need cyclic pi-conjugation
  R5  < 6 heavy atoms                          -> too small (O2, etc.)
  R6  < 6 sp2/aromatic (conjugated) atoms      -> insufficient pi-conjugation

Nothing here fabricates data: it only *removes* rows and re-ranks the rest.
Run:  python filter_genuine_osc.py
"""
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors  # noqa: F401  (kept for interactive use)
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

CSV = "../../data/PCE_paper_GDB9.csv"

# Metals / metalloids that mark a row as inorganic / organometallic salt.
METALS = {
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "Rb", "Sr", "Y", "Zr", "Nb",
    "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Cs", "Ba",
    "La", "Ce", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "B", "Si", "As", "Te",
}


def reject_reason(smiles):
    """Return None if the molecule is a plausible OSC, else a short reason."""
    if not isinstance(smiles, str) or not smiles:
        return "unparsable"
    if "." in smiles:
        return "R1_multifragment"
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "unparsable"
    syms = {a.GetSymbol() for a in mol.GetAtoms()}
    if syms & METALS:
        return "R2_metal"
    ri = mol.GetRingInfo()
    if ri.NumRings() == 0:
        return "R3_no_ring"
    if not any(a.GetIsAromatic() for a in mol.GetAtoms()):
        return "R4_no_aromatic_ring"
    heavy = mol.GetNumHeavyAtoms()
    if heavy < 6:
        return "R5_too_small"
    conj = sum(
        1
        for a in mol.GetAtoms()
        if a.GetIsAromatic() or a.GetHybridization() == Chem.HybridizationType.SP2
    )
    if conj < 6:
        return "R6_low_conjugation"
    return None


def main():
    df = pd.read_csv(CSV)
    df["PCE_max"] = df[["pce_pcbm(%)", "pce_pcdtbt(%)"]].max(axis=1)
    df["PCE_SAScore"] = df["PCE_max"] - df["sas1(%)"]
    df["reject"] = df["SMILES"].map(reject_reason)

    kept = df[df["reject"].isna()].sort_values("PCE_SAScore", ascending=False)
    rejected = df[df["reject"].notna()]

    print(f"Total molecules           : {len(df)}")
    print(f"Kept (genuine OSC)        : {len(kept)}")
    print(f"Rejected (artifacts)      : {len(rejected)}")
    print()
    print("Rejection breakdown:")
    print(rejected["reject"].value_counts().to_string())
    print()

    print("=== REAL TOP-15 genuine organic-semiconductor candidates ===")
    cols = ["mol_id", "formula", "PCE_SAScore", "sas1(%)", "PCE_max", "SMILES"]
    top = kept.head(15)[cols]
    for _, r in top.iterrows():
        print(
            f"  {int(r['mol_id']):>6}  PS={r['PCE_SAScore']:6.2f}  "
            f"SAS={r['sas1(%)']:.2f}  PCE={r['PCE_max']:5.1f}  "
            f"{r['formula']:<14} {r['SMILES']}"
        )

    print()
    print("=== Fate of the OLD (contaminated) top-7 ===")
    for mid in [17851, 20778, 4550, 977, 11029, 1712, 7801]:
        row = df[df["mol_id"] == mid].iloc[0]
        verdict = "KEPT" if pd.isna(row["reject"]) else f"REJECTED ({row['reject']})"
        print(f"  {mid:>6} {row['formula']:<12} PS={row['PCE_SAScore']:6.2f}  -> {verdict}")

    out = "genuine_osc_ranked.csv"
    kept[cols].to_csv(out, index=False)
    print(f"\nWrote full ranked survivor list -> {out}")


if __name__ == "__main__":
    sys.exit(main())
