"""
B3LYP/6-31G* geometry optimization for the 3 NEW gap-filling candidates
(20549 stable, 5081 reactive/azide, 23424 reactive/nitroso), added 2026-08-03
to break the heavy-atom/class-separation confound flagged by the JCIM desk
rejection (stable donors 26-33 heavy atoms vs reactive structures 15-19, no
overlap). 20549 and 5081 both have 24 heavy atoms; 23424 has 23 -- all inside
the former empty gap.

Same pipeline as `optimize_new_geometries.py` (used for 17574/18506/9168 on
2026-07-16): RDKit ETKDGv3(seed=42)+UFF conformer -> ORCA
"! B3LYP 6-31G* Opt TightSCF" -> optimized .xyz into geometries/, so all 10
candidate molecules (original 7 + these 3) share one geometric footing.

Usage (from this directory, venv active):
    python3 optimize_gap_geometries.py
    python3 optimize_gap_geometries.py --parse-only   # re-parse existing .out
"""
import argparse
import os
import subprocess
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

HOME = Path.home()
ORCA_EXE = HOME / "orca-6-1-1" / "orca"

BASE = Path(__file__).resolve().parent
GEOM_DIR = BASE / "geometries"
OPT_DIR = BASE / "opt_gap_geometries"

SEED = 42
CHARGE = 0
MULTIPLICITY = 1


def _physical_cores() -> int:
    try:
        core_ids = set()
        with open("/proc/cpuinfo") as fh:
            phys = core = None
            for line in fh:
                if line.startswith("physical id"):
                    phys = line.split(":")[1].strip()
                elif line.startswith("core id"):
                    core = line.split(":")[1].strip()
                elif not line.strip() and phys is not None and core is not None:
                    core_ids.add((phys, core))
                    phys = core = None
        if core_ids:
            return len(core_ids)
    except OSError:
        pass
    return max(1, (os.cpu_count() or 1) // 2)


DEFAULT_NPROCS = _physical_cores()
MAXCORE_MB = 3500

# From genuine_osc_ranked.csv (chemical-validity filter applied, same as the
# original 7); SMILES verbatim from that CSV.
MOLECULES = {
    "20549": "OCCNc1ccc(c2c1C(=O)c1c(C2=O)cccc1)NCCO",
    "5081": "CCOC(=O)c1ncn2c1CN(C)C(=O)c1c2ccc(c1)N=[N+]=[N-]",
    "23424": "O=Nc1ccc(cc1)N(CCOS(=O)(=O)C)CCOS(=O)(=O)C",
}


def rdkit_conformer_xyz(mol_id: str, smiles: str, run_subdir: Path) -> Path:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit failed to parse SMILES for {mol_id}: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError(f"RDKit embedding failed for {mol_id}")
    AllChem.UFFOptimizeMolecule(mol, maxIters=2000)

    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), f"{mol_id} RDKit ETKDGv3(seed={SEED})+UFF"]
    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol():<2} {p.x:>14.8f} {p.y:>14.8f} {p.z:>14.8f}")
    xyz_path = run_subdir / f"{mol_id}_guess.xyz"
    xyz_path.write_text("\n".join(lines) + "\n")
    return xyz_path


def write_opt_input(mol_id: str, guess_xyz: Path, nprocs: int, run_subdir: Path) -> Path:
    inp_path = run_subdir / f"{mol_id}_opt.inp"
    inp_path.write_text(
        f"! B3LYP 6-31G* Opt TightSCF\n"
        f"%pal nprocs {nprocs} end\n"
        f"%maxcore {MAXCORE_MB}\n"
        f"* xyzfile {CHARGE} {MULTIPLICITY} {guess_xyz.name}\n"
    )
    return inp_path


def run_orca(inp_path: Path) -> Path:
    out_path = inp_path.with_suffix(".out")
    with open(out_path, "w") as fh:
        subprocess.run(
            [str(ORCA_EXE), str(inp_path)],
            stdout=fh, stderr=subprocess.STDOUT, cwd=inp_path.parent, check=False,
        )
    return out_path


def parse_optimized_xyz(mol_id: str, run_subdir: Path, out_path: Path) -> Path:
    text = out_path.read_text(errors="replace")
    if "ORCA TERMINATED NORMALLY" not in text:
        raise RuntimeError(f"ORCA did not terminate normally for {mol_id}: {out_path}")
    if "THE OPTIMIZATION HAS CONVERGED" not in text:
        raise RuntimeError(f"Optimization did not converge for {mol_id}: {out_path}")
    orca_xyz = run_subdir / f"{mol_id}_opt.xyz"
    if not orca_xyz.exists():
        raise FileNotFoundError(f"ORCA optimized xyz missing for {mol_id}: {orca_xyz}")
    dest = GEOM_DIR / f"{mol_id}.xyz"
    dest.write_text(orca_xyz.read_text())
    return dest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--parse-only", action="store_true")
    parser.add_argument("--nprocs", type=int, default=DEFAULT_NPROCS)
    args = parser.parse_args()

    OPT_DIR.mkdir(exist_ok=True)
    GEOM_DIR.mkdir(exist_ok=True)

    for mol_id, smiles in MOLECULES.items():
        run_subdir = OPT_DIR / mol_id
        run_subdir.mkdir(parents=True, exist_ok=True)
        inp_path = run_subdir / f"{mol_id}_opt.inp"
        out_path = inp_path.with_suffix(".out")

        if not args.parse_only:
            if not ORCA_EXE.exists():
                raise FileNotFoundError(
                    f"ORCA not found at {ORCA_EXE}. Paths resolve from "
                    f"Path.home() ({HOME}) -- check orca-6-1-1 exists there."
                )
            guess_xyz = rdkit_conformer_xyz(mol_id, smiles, run_subdir)
            write_opt_input(mol_id, guess_xyz, args.nprocs, run_subdir)
            print(f"Optimizing {mol_id} (B3LYP/6-31G*) ...")
            run_orca(inp_path)

        dest = parse_optimized_xyz(mol_id, run_subdir, out_path)
        print(f"  {mol_id}: optimized geometry written to {dest}")

    print("\nAll three gap-filling geometries optimized and placed in geometries/.")


if __name__ == "__main__":
    main()
