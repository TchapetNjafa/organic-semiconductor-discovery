"""sTDA-xTB S1 for HOPV15 with the exact shard protocol, for calibration against experimental optical gaps.

Protocol (identical to NTO_and_related_study/scripts/00 + 01): SMILES -> RDKit AddHs -> ETKDG(seed 42)
-> UFF(maxIters=200) -> `xtb4stda mol.xyz` single point -> `stda -xtb wfn.xtb` (defaults).

usage: python -I hopv_stda.py <HOPV_15_revised_2.data> <out_dir> [n_jobs]
writes <out_dir>/hopv15_records.csv and <out_dir>/hopv15_stda.csv
"""
import csv
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

HARTREE_EV = 27.211386245988
EXP_COLS = ["doi", "inchikey", "construction", "architecture", "complement", "homo_exp", "lumo_exp",
            "gap_echem", "gap_opt", "pce", "voc", "jsc", "ff"]
XTB4STDA = os.path.expanduser("~/xtb4stda/xtb4stda")
STDA = os.path.expanduser("~/stda_v1-6/stda")


def parse_hopv(path):
    """One dict per molecule: experimental line + B3LYP/def2-SVP gap of conformer 1."""
    lines = Path(path).read_text().splitlines()
    out, i = [], 0
    while i < len(lines):
        smiles, exp = lines[i].strip(), lines[i + 2].split(",")
        rec = {"smiles": smiles, **dict(zip(EXP_COLS, (x.strip() for x in exp)))}
        n_conf = int(lines[i + 4])
        i += 5
        for c in range(n_conf):
            n_at = int(lines[i + 1])
            dft = lines[i + 2 + n_at: i + 6 + n_at]
            if c == 0:
                b3 = next(d for d in dft if d.startswith("QChem B3LYP")).split(",")
                rec["homo_b3lyp_eV"] = float(b3[1]) * HARTREE_EV
                rec["lumo_b3lyp_eV"] = float(b3[2]) * HARTREE_EV
                rec["gap_b3lyp_eV"] = float(b3[3]) * HARTREE_EV
            i += 2 + n_at + 4
        out.append(rec)
    return out


def stda_s1(job):
    idx, smiles, out_dir = job
    wd = Path(out_dir) / f"hopv_{idx:03d}"
    wd.mkdir(parents=True, exist_ok=True)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return idx, "invalid_smiles", None, None, None, None
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        return idx, "embed_failed", None, None, None, None
    AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    Chem.MolToXYZFile(mol, str(wd / "mol.xyz"))
    env = {**os.environ, "XTB4STDAHOME": os.path.expanduser("~/xtb4stda"), "OMP_NUM_THREADS": "1"}
    for stale in ("wfn.xtb", "tda.dat"):
        (wd / stale).unlink(missing_ok=True)
    try:
        a =subprocess.run([XTB4STDA, "mol.xyz"], cwd=wd, env=env, capture_output=True, timeout=1800)
        b = subprocess.run([STDA, "-xtb", "wfn.xtb"], cwd=wd, env=env, capture_output=True, timeout=1800)
        if a.returncode or b.returncode:
            return (idx, "nonzero_exit", None, None, None, None)
    except subprocess.TimeoutExpired:
        return idx, "timeout", None, None, None, None
    tda = wd / "tda.dat"
    if not tda.exists():
        return idx, "stda_failed", None, None, None, None
    (wd / "wfn.xtb").unlink(missing_ok=True)
    rows = [l.split() for l in tda.read_text().split("DATXY", 1)[1].splitlines() if len(l.split()) >= 3]
    s1 = rows[0] if rows else None
    br = next((r for r in rows if float(r[2]) > 0.01), None)
    if s1 is None:
        return idx, "no_states", None, None, None, None
    return (idx, "success", float(s1[1]), float(s1[2]),
            float(br[1]) if br else None, float(br[2]) if br else None)


def main():
    src, out_dir = sys.argv[1], sys.argv[2]
    n_jobs = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    recs = parse_hopv(src)
    assert len(recs) == 350 and recs[0]["gap_opt"] == "1.59", "HOPV15 parse drifted from the v4 file"
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    with open(Path(out_dir) / "hopv15_records.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["idx", *recs[0].keys()])
        w.writeheader()
        for k, r in enumerate(recs):
            w.writerow({"idx": k, **r})
    jobs = [(k, r["smiles"], out_dir) for k, r in enumerate(recs)]
    with ProcessPoolExecutor(n_jobs) as ex, open(Path(out_dir) / "hopv15_stda.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["idx", "status", "E_S1_eV", "f_S1", "E_bright_eV", "f_bright"])
        for row in ex.map(stda_s1, jobs):
            w.writerow(row)
            f.flush()


if __name__ == "__main__":
    main()
