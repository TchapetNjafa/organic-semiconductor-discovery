#!/usr/bin/env python3
"""
Reactive-group screen on top of the chemical-validity filter (rule R7).

Reads the chemically valid molecules written by filter_genuine_osc.py and flags
reactive or energetic functional groups that make a molecule an unrealistic OPV
material even though it is a valid conjugated structure:

  azide    [#7]-[#7+]#[#7] / [#7]=[#7+]=[#7-] / [#7-]-[#7+]#[#7]  -> energetic
  nitroso  [#6][NX2]=[OX1]                                        -> reactive
  diazo    [#6]=[#7+]=[#7-]                                       -> reactive

This is an inspection diagnostic, NOT a safety filter. It covers three motifs.
It does not test for peroxides, nitro groups, perchlorates, strained rings,
N-oxides, gem-polynitro or acyl azides, and a molecule that passes it may still
be hazardous or photochemically unstable. The main text's limitations section
says so explicitly; anyone reusing this file should read that before treating a
clean result as a stability claim.

Two views of the region above threshold are produced:
  View A  chemical validity only        (reactive groups retained)
  View B  chemical validity + this screen (reactive groups removed)

No data is fabricated: molecules are only flagged and, in View B, removed.

Run after filter_genuine_osc.py:  python stability_screen.py
"""
import argparse
import pathlib
import sys

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = pathlib.Path(__file__).resolve().parent

# Each pattern is written atom by atom so it can be lifted into another screen
# without re-deriving it. [#7] is any nitrogen, [#6] any carbon; a charge in the
# bracket is a formal charge, and NX2/OX1 fix the connectivity (two connections
# on N, one on O) so the pattern cannot match a larger group that merely
# contains the same atoms.
#
# The three azide patterns are one motif written in the three protonation- and
# resonance-forms that databases actually store; a molecule is flagged if it
# matches any of them. This is the same tautomer/resonance bookkeeping problem
# the main text raises for o-quinones: a SMARTS matches the form a structure was
# written in, not the compound. Canonicalising with
# rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumerator before matching
# is the general fix; here we enumerate the forms by hand instead.
REACTIVE = {
    # Organic azide, R-N3. Energetic: the terminal N2 leaves under heat, light
    # or shock, which rules the molecule out as a photovoltaic absorber.
    "azide": [
        "[#7]-[#7+]#[#7]",     # neutral-N / N+ triple-bonded N: R-N=[N+]=[N-] drawn linear
        "[#7]=[#7+]=[#7-]",    # the cumulated all-double-bond resonance form
        "[#7-]-[#7+]#[#7]",    # anionic form, e.g. an azide written as a separate ion
    ],
    # C-nitroso, R-N=O. Reactive: dimerises, oxidises, and absorbs in the
    # visible, so it is neither photostable nor spectrally innocent.
    # NX2 = nitrogen with exactly two connections, OX1 = terminal oxygen, which
    # together exclude nitro (R-NO2) and N-oxides from this pattern.
    "nitroso": ["[#6][NX2]=[OX1]"],
    # Diazo, R2C=[N+]=[N-]. Loses N2 to give a carbene; used as a reagent
    # precisely because it is unstable.
    "diazo": ["[#6]=[#7+]=[#7-]"],
}
PATTS = {name: [Chem.MolFromSmarts(s) for s in sm] for name, sm in REACTIVE.items()}

# Motifs this screen does NOT cover, named so the omission is on the record.
NOT_SCREENED = [
    "peroxides and hydroperoxides", "nitro and gem-polynitro groups",
    "perchlorates and other oxidising counter-ions", "strained rings",
    "N-oxides", "acyl azides and azidoformates", "tetrazoles and triazoles",
    "diazonium salts", "fulminates and nitrile oxides",
]


def flags(smiles):
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None:
        return []
    return [name for name, patts in PATTS.items()
            if any(mol.HasSubstructMatch(p) for p in patts)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--input", type=pathlib.Path,
                    default=HERE / "genuine_osc_ranked.csv")
    ap.add_argument("--score", default="PCE_SAScore_compl_FF065")
    a = ap.parse_args()

    if not a.input.exists():
        sys.exit(f"missing {a.input}\nRun filter_genuine_osc.py first.")

    df = pd.read_csv(a.input)
    df["reactive_recomputed"] = df["SMILES"].map(lambda s: ",".join(flags(s)))

    # Cross-check against the flags scharber_recompute.py wrote.
    pipe = df["reactive_groups"].fillna("").astype(str)
    bad = df.loc[pipe != df["reactive_recomputed"],
                 ["mol_id", "reactive_groups", "reactive_recomputed"]]
    if len(bad):
        print("PIPELINE DISAGREEMENT on", len(bad), "rows:")
        print(bad.head(10).to_string(index=False))
        sys.exit("reactive-group SMARTS do not match scharber_recompute.py")
    print(f"SMARTS cross-check: {len(df)} rows agree with scharber_recompute.py")

    df["stable_recomputed"] = df["reactive_recomputed"] == ""
    admitted = df[df[a.score].notna()].sort_values(a.score, ascending=False)
    viable = admitted[admitted[a.score] > 0]
    view_a = viable
    view_b = viable[viable["stable_recomputed"]]

    def show(title, sub):
        print(f"\n=== {title}  ({len(sub)} molecules) ===")
        for _, r in sub.iterrows():
            rg = r["reactive_recomputed"] or "-"
            print(f"  {int(r['mol_id']):>6}  PS={r[a.score]:7.2f}  "
                  f"SAS={r['SAScore']:.2f}  Eg={r['gap_eV']:.3f}  "
                  f"[{rg:<8}] {r['formula']:<14} {r['SMILES']}")

    print(f"\nChemically valid molecules the gates admit: {len(admitted)}")
    print(f"  of which above threshold ({a.score} > 0): {len(viable)}")
    show(f"VIEW A  validity only, reactive groups RETAINED", view_a)
    show(f"VIEW B  validity + reactive-group screen", view_b)

    print("\n=== Reactive-group flags among the molecules above threshold ===")
    for _, r in view_a.iterrows():
        print(f"  {int(r['mol_id']):>6} {r['formula']:<14} -> "
              f"{r['reactive_recomputed'] or 'none flagged'}")

    print("\n=== Motifs this screen does NOT test for ===")
    for m in NOT_SCREENED:
        print(f"  - {m}")
    print("A clean result here is not a stability or safety claim.")

    view_a.to_csv(HERE / "viable_view_a_metric.csv", index=False)
    view_b.to_csv(HERE / "viable_view_b_stability_screened.csv", index=False)
    print("\nWrote viable_view_a_metric.csv and viable_view_b_stability_screened.csv")


if __name__ == "__main__":
    sys.exit(main())
