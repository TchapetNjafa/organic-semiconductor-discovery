#!/usr/bin/env python3
"""
Stability screen on top of the chemical-validity filter.

Reads the genuine-OSC survivors (genuine_osc_ranked.csv) and flags
reactive / unstable functional groups that make a molecule an unrealistic
real-world OPV material even though it is a valid conjugated structure:

  azide   [#7]-[#7+]#[#7]  or  [#7]=[#7+]=[#7-]   -> energetic / explosophoric
  nitroso [#6][NX2]=[OX1]                          -> reactive, often unstable
  diazo   [#6]=[#7+]=[#7-]                          -> reactive

Produces two honest "views" of the viable region (PCE_SAScore > 0):
  View A  chemical-validity filter only          (azides included)
  View B  chemical-validity + stability screen    (reactive groups removed)

No data is fabricated: molecules are only flagged and, in View B, removed.
Run after filter_genuine_osc.py:  python stability_screen.py
"""
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

REACTIVE = {
    "azide": ["[#7]-[#7+]#[#7]", "[#7]=[#7+]=[#7-]", "[#7-]-[#7+]#[#7]"],
    "nitroso": ["[#6][NX2]=[OX1]"],
    "diazo": ["[#6]=[#7+]=[#7-]"],
}
PATTS = {name: [Chem.MolFromSmarts(s) for s in sm] for name, sm in REACTIVE.items()}


def flags(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    hit = []
    for name, patts in PATTS.items():
        if any(mol.HasSubstructMatch(p) for p in patts):
            hit.append(name)
    return hit


def main():
    df = pd.read_csv("genuine_osc_ranked.csv")
    df["reactive_groups"] = df["SMILES"].map(lambda s: ",".join(flags(s)))
    df["stable"] = df["reactive_groups"] == ""

    viable = df[df["PCE_SAScore"] > 0].copy()
    view_a = viable
    view_b = viable[viable["stable"]]

    def show(title, sub):
        print(f"\n=== {title}  ({len(sub)} molecules, PCE_SAScore > 0) ===")
        for _, r in sub.iterrows():
            rg = r["reactive_groups"] or "-"
            print(
                f"  {int(r['mol_id']):>6}  PS={r['PCE_SAScore']:6.2f}  "
                f"SAS={r['sas1(%)']:.2f}  {r['formula']:<14} "
                f"[{rg:<8}] {r['SMILES']}"
            )

    show("VIEW A  chemical-validity filter (azides INCLUDED)", view_a)
    show("VIEW B  + stability screen (reactive groups REMOVED)", view_b)

    print("\n=== Reactive-group flags among the 4 metric-viable molecules ===")
    for _, r in view_a.iterrows():
        print(f"  {int(r['mol_id']):>6} {r['formula']:<14} -> "
              f"{r['reactive_groups'] or 'none (stable)'}")

    view_a.to_csv("viable_view_a_metric.csv", index=False)
    view_b.to_csv("viable_view_b_stability_screened.csv", index=False)
    print("\nWrote viable_view_a_metric.csv and viable_view_b_stability_screened.csv")


if __name__ == "__main__":
    main()
