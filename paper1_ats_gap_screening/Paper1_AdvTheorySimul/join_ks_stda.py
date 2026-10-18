"""Join the stored PubChemQC Kohn-Sham gap with the sTDA-xTB energies (inner join on mol_id).

usage: python -I join_ks_stda.py <project_root> [stda_csv] [out_csv]   -> outputs/extract/ks_vs_stda_joined.csv
delta = E_S1 - E_KS [eV].
"""
import sys
from pathlib import Path

import pandas as pd

if __name__ == "__main__":
    root = Path(sys.argv[1])
    stda = Path(sys.argv[2]) if len(sys.argv) > 2 else root / "outputs/extract/stda_all.csv"
    dest = Path(sys.argv[3]) if len(sys.argv) > 3 else root / "outputs/extract/ks_vs_stda_joined.csv"
    d = pd.read_csv(root / "ZENODO_UPLOAD/dataset_pubchemqc_opv_17458.csv", usecols=["mol_id", "formula", "GAP(eV)"])
    s = pd.read_csv(stda, dtype=str).apply(pd.to_numeric, errors="coerce").dropna(subset=["mol_id"])
    s["mol_id"] = s.mol_id.astype(int)
    j = d.merge(s, on="mol_id", how="inner")
    j["delta"] = j.E_S1_eV - j["GAP(eV)"]
    j.to_csv(dest, index=False)
    print(len(j), "rows ->", dest)
