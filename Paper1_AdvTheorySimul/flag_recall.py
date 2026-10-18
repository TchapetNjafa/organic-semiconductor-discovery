"""Flag recall: B3LYP/6-31G* recompute on a random sample of UNflagged single-fragment records (L3, death-angle A1/A5).

usage: python -I flag_recall.py prepare <root> [n=120] | run <root> <shard> <nshards>
Defect = recomputed gap exceeds the stored gap by > 0.5 eV (same criterion and UFF geometry protocol as L-2026-10-17-09).
Resumable: results appended to outputs/extract/flag-recall/result_<shard>.tsv; finished mol_ids are skipped.
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

root = Path(sys.argv[2])
out = root / "outputs/extract/flag-recall"
HARTREE_EV = 27.211386245988


def prepare(n):
    from rdkit import Chem
    j = pd.read_csv(root / "outputs/extract/ks_vs_stda_joined.csv")
    smi = pd.read_csv(root / "outputs/extract/gate_cascade.csv", usecols=["mol_id", "SMILES"])
    j = j.merge(smi, on="mol_id")
    used = set(pd.read_csv(root / "outputs/extract/record-validity-result.csv").mol_id)
    j = j[~j.SMILES.str.contains(r"\.") & ((j.E_S1_eV - j["GAP(eV)"]) <= 1.0) & ~j.mol_id.isin(used)]
    j = j[j.SMILES.map(lambda s: (m := Chem.MolFromSmiles(s)) is not None and 6 <= m.GetNumHeavyAtoms() <= 35)]
    s = j.sample(n=n, random_state=11).reset_index(drop=True)
    s.to_csv(out / "sample.csv", index=False)
    print("universe", len(j), "sampled", len(s))


def xyz(smiles, path):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(m, randomSeed=42) != 0:
        return False
    AllChem.UFFOptimizeMolecule(m, maxIters=500)
    Chem.MolToXYZFile(m, str(path))
    return True


def run(shard, nsh):
    from pyscf import dft, gto, lib
    lib.num_threads(int(os.environ.get("NTHREADS", "4")))
    s = pd.read_csv(out / "sample.csv")
    resf = out / f"result_{shard}.tsv"
    done = set(pd.read_csv(resf, sep="\t").mol_id) if resf.exists() else set()
    for i, r in s.iloc[shard::nsh].iterrows():
        if r.mol_id in done:
            continue
        p = out / "xyz" / f"{r.mol_id}.xyz"
        if not xyz(r.SMILES, p):
            line = f"{r.mol_id}\tembed_fail\tnan\tnan\n"
        else:
            mf = dft.RKS(gto.M(atom=str(p), basis="6-31g*", verbose=0), xc="b3lyp").run()
            e, o = mf.mo_energy, mf.mo_occ
            g = (e[o == 0].min() - e[o > 0].max()) * HARTREE_EV
            line = f"{r.mol_id}\t{'ok' if mf.converged else 'noconv'}\t{r['GAP(eV)']:.4f}\t{g:.4f}\n"
        new = not resf.exists()
        with open(resf, "a") as f:
            if new:
                f.write("mol_id\tstatus\tgap_stored\tgap_recomp\n")
            f.write(line)


if sys.argv[1] == "prepare":
    prepare(int(sys.argv[3]) if len(sys.argv) > 3 else 120)
else:
    run(int(sys.argv[3]), int(sys.argv[4]))
