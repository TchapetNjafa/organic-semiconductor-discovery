"""
Dock three ADDITIONAL stable, validity-passing donor candidates
(17574, 18506, 9168) against the same four receptors and under the identical
protocol used for the original four candidates in redock_candidates.py.

These three molecules pass the same chemical-validity + reactive-group
stability screen as 1712 but were below the PCE_SAScore > 0 cutoff in the
companion photovoltaic screen. They are added here to broaden the docking set
beyond a single stable candidate and to populate the middle of the
dipole-versus-binding trade-off axis with genuinely stable structures.

Protocol is byte-for-byte the same as redock_candidates.py:
- smina by full path; RDKit ETKDGv3 (seed 42) + UFF; single obabel PDBQT pass
  (gasteiger charges, no second 3D embedding); whole-protein blind-docking box
  from receptor heavy-atom extent + 8 A padding; fixed exhaustiveness/seed/modes.

Appends its rows to docking_results_clean.csv (does not overwrite the
original four candidates), and writes a standalone
docking_results_new_candidates.csv for provenance.
"""
import csv
import subprocess
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

SMINA = "/home/tchapet/SMINA/smina.static"
BASE = Path(__file__).resolve().parent
RECEPTOR_DIR = Path(
    "/home/tchapet/Post-Doc/ARTICLES-DOSSIERS-DES-ARTICLES-EN-REDACTION/"
    "NOUVELS-AXES-DE-RECHERCHE-A-REGARDER-URGEMMENT/ARTICLES-EN-REDACTIONS/"
    "ARTICLE_MVOTO/MVOTOOPV/PUBCHEMQC_MULTIFONCTIONAL_2025/receptors"
)
LIGAND_DIR = BASE / "ligands"
OUT_DIR = BASE / "outputs"
LOG_DIR = BASE / "logs"

SEED = 42
EXHAUSTIVENESS = 32
NUM_MODES = 9
PADDING_ANGSTROM = 8.0

# New stable, validity-passing donors (SMILES from genuine_osc_ranked.csv).
LIGANDS = {
    "17574": "CC(=O)N(c1ccc(cc1)Nc1cc(c(c2c1C(=O)c1ccccc1C2=O)N)S(=O)(=O)O)C",
    "18506": "CCCCn1c(=O)c2c(c1=O)c(N)c1c(c2N)c(=O)c2c(c1=O)cccc2",
    "9168": "c1ccc2c(c1)c1ccc3c(c1cc2)ccc1c3ccc2c1cccc2",
}

RECEPTORS = {
    "HIV1protease": "1DMP",
    "Hsp90": "2XJX",
    "Neurodegenerative_1SYH": "1SYH",
    "COVID19_6Y2F": "6Y2F",
}


def prepare_ligand(mol_id: str, smiles: str) -> Path:
    pdb_path = LIGAND_DIR / f"{mol_id}.pdb"
    pdbqt_path = LIGAND_DIR / f"{mol_id}.pdbqt"

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit failed to parse SMILES for {mol_id}: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    embed_status = AllChem.EmbedMolecule(mol, params)
    if embed_status != 0:
        raise RuntimeError(f"RDKit embedding failed for {mol_id}")
    AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
    Chem.MolToPDBFile(mol, str(pdb_path))

    subprocess.run(
        [
            "obabel",
            str(pdb_path),
            "-O", str(pdbqt_path),
            "--partialcharge", "gasteiger",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return pdbqt_path


def receptor_box(receptor_pdb: Path) -> dict:
    xs, ys, zs = [], [], []
    with open(receptor_pdb) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    xs.append(float(line[30:38]))
                    ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))
                except ValueError:
                    continue
    center = {
        "x": (max(xs) + min(xs)) / 2,
        "y": (max(ys) + min(ys)) / 2,
        "z": (max(zs) + min(zs)) / 2,
    }
    size = {
        "x": (max(xs) - min(xs)) + 2 * PADDING_ANGSTROM,
        "y": (max(ys) - min(ys)) + 2 * PADDING_ANGSTROM,
        "z": (max(zs) - min(zs)) + 2 * PADDING_ANGSTROM,
    }
    return {"center": center, "size": size}


def run_smina(ligand_pdbqt: Path, receptor_pdbqt: Path, box: dict, out_pdbqt: Path, log_path: Path) -> float:
    cmd = [
        SMINA,
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", str(box["center"]["x"]),
        "--center_y", str(box["center"]["y"]),
        "--center_z", str(box["center"]["z"]),
        "--size_x", str(box["size"]["x"]),
        "--size_y", str(box["size"]["y"]),
        "--size_z", str(box["size"]["z"]),
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", str(NUM_MODES),
        "--seed", str(SEED),
        "--out", str(out_pdbqt),
        "--log", str(log_path),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"smina failed: {result.stderr}\n{result.stdout}")

    with open(log_path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 2 and parts[0] == "1":
                try:
                    return float(parts[1])
                except ValueError:
                    continue
    raise RuntimeError(f"could not parse best-mode affinity from {log_path}")


def main():
    ligand_pdbqts = {mid: prepare_ligand(mid, smi) for mid, smi in LIGANDS.items()}

    boxes = {}
    for target, pdb_id in RECEPTORS.items():
        boxes[target] = receptor_box(RECEPTOR_DIR / f"{pdb_id}.pdb")

    rows = []
    for target, pdb_id in RECEPTORS.items():
        receptor_pdbqt = RECEPTOR_DIR / f"{pdb_id}.pdbqt"
        box = boxes[target]
        for mol_id, ligand_pdbqt in ligand_pdbqts.items():
            out_pdbqt = OUT_DIR / f"{mol_id}_{target}.pdbqt"
            log_path = LOG_DIR / f"{mol_id}_{target}.log"
            affinity = run_smina(ligand_pdbqt, receptor_pdbqt, box, out_pdbqt, log_path)
            print(f"{mol_id} vs {target} ({pdb_id}): {affinity} kcal/mol")
            rows.append({
                "mol_id": mol_id,
                "target": target,
                "pdb_id": pdb_id,
                "affinity_kcal_mol": affinity,
                "box_center_x": box["center"]["x"],
                "box_center_y": box["center"]["y"],
                "box_center_z": box["center"]["z"],
                "box_size_x": box["size"]["x"],
                "box_size_y": box["size"]["y"],
                "box_size_z": box["size"]["z"],
                "exhaustiveness": EXHAUSTIVENESS,
                "seed": SEED,
                "num_modes": NUM_MODES,
            })

    # Standalone provenance file.
    new_path = BASE / "docking_results_new_candidates.csv"
    with open(new_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {len(rows)} results to {new_path}")

    # Append to the combined clean file (idempotent: rebuild from scratch to
    # avoid duplicate rows if re-run).
    combined_path = BASE / "docking_results_clean.csv"
    original = []
    with open(combined_path) as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames
        for r in reader:
            if r["mol_id"] not in LIGANDS:  # drop any prior copies of the new set
                original.append(r)
    with open(combined_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(original)
        writer.writerows(rows)
    print(f"Appended {len(rows)} rows to {combined_path} "
          f"({len(original)} original + {len(rows)} new)")


if __name__ == "__main__":
    main()
