"""
Apo-receptor docking of the 3 heavy-atom-gap-filling candidates (JCIM
desk-rejection fix #4, 2026-08-03): 20549 (stable, C18H18N2O4, 24 heavy
atoms), 5081 (reactive/azide, C15H14N6O3, 24 heavy atoms -- exact size
match to 20549), 23424 (reactive/nitroso, C12H18N2O7S2, 23 heavy atoms).

These were selected from the 10,215-molecule genuine_osc_ranked.csv pool
(validity-filtered PubChemQC subset) specifically to break the perfect
size/class confound the JCIM associate editor flagged (original 7-candidate
set: stable = 26/27/31/33 heavy atoms, reactive = 15/16/19 heavy atoms, a
total unoccupied gap at 20-25). See analysis-ledger.md entry
L-2026-08-03-02 for the selection quantification.

Protocol is identical to dock_new_candidates.py / redock_candidates_apo.py:
RDKit ETKDGv3 (seed=42) + UFF conformer -> PDB -> PDBQT via obabel
(gasteiger charges, single embedding, no --gen3d) -> whole-protein blind
docking against the SAME apo receptors used for the apo-redocking fix
(calibration_outputs/{pdb_id}_apo.pdbqt), same box (from calibration csv),
exhaustiveness=32, num_modes=9, seed=42.

Appends rows to docking_results_apo.csv (does not duplicate/overwrite the
7 candidates already redocked there) and writes a standalone
docking_results_gap_candidates_apo.csv for provenance.
"""
import csv
import subprocess
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

SMINA = "/home/tchapet/SMINA/smina.static"
BASE = Path(__file__).resolve().parent
LIGAND_DIR = BASE / "ligands"
APO_DIR = BASE / "calibration_outputs"
OUT_DIR = BASE / "outputs_apo"
LOG_DIR = BASE / "logs_apo"
LIGAND_DIR.mkdir(exist_ok=True)
OUT_DIR.mkdir(exist_ok=True)
LOG_DIR.mkdir(exist_ok=True)

SEED = 42
EXHAUSTIVENESS = 32
NUM_MODES = 9

# SMILES from genuine_osc_ranked.csv (Paper1_Digital_Discovery/), verified
# by RDKit heavy-atom count against the ledger entry before use.
LIGANDS = {
    "20549": "OCCNc1ccc(c2c1C(=O)c1c(C2=O)cccc1)NCCO",
    "5081": "CCOC(=O)c1ncn2c1CN(C)C(=O)c1c2ccc(c1)N=[N+]=[N-]",
    "23424": "O=Nc1ccc(cc1)N(CCOS(=O)(=O)C)CCOS(=O)(=O)C",
}

# Same 4 targets, same apo receptors/box already used in redock_candidates_apo.py.
TARGETS = {
    "HIV1protease": ("1DMP", {
        "center": {"x": -11.782, "y": 20.5045, "z": 26.767999999999997},
        "size": {"x": 59.028, "y": 54.29900000000001, "z": 69.70400000000001},
    }),
    "Hsp90": ("2XJX", {
        "center": {"x": 12.0865, "y": -2.5024999999999995, "z": 21.817},
        "size": {"x": 67.645, "y": 63.191, "z": 71.034},
    }),
    "Neurodegenerative_1SYH": ("1SYH", {
        "center": {"x": 17.6955, "y": 11.253, "z": 23.977500000000003},
        "size": {"x": 70.87299999999999, "y": 73.434, "z": 73.31700000000001},
    }),
    "COVID19_6Y2F": ("6Y2F", {
        "center": {"x": -5.699999999999999, "y": -0.010999999999999233, "z": 13.115},
        "size": {"x": 69.406, "y": 82.374, "z": 71.038},
    }),
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

    rows = []
    for target, (pdb_id, box) in TARGETS.items():
        apo_receptor = APO_DIR / f"{pdb_id}_apo.pdbqt"
        if not apo_receptor.exists():
            raise FileNotFoundError(f"missing apo receptor: {apo_receptor}")
        for mol_id, ligand_pdbqt in ligand_pdbqts.items():
            out_pdbqt = OUT_DIR / f"{mol_id}_{target}_apo.pdbqt"
            log_path = LOG_DIR / f"{mol_id}_{target}_apo.log"
            affinity = run_smina(ligand_pdbqt, apo_receptor, box, out_pdbqt, log_path)
            print(f"{mol_id} vs {target} (apo {pdb_id}): {affinity} kcal/mol")
            rows.append({
                "mol_id": mol_id,
                "target": target,
                "pdb_id": pdb_id,
                "receptor_state": "apo",
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

    standalone_path = BASE / "docking_results_gap_candidates_apo.csv"
    with open(standalone_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {len(rows)} results to {standalone_path}")

    # Append to the combined apo file (idempotent: rebuild, dropping any
    # prior copies of these 3 mol_ids, so re-runs don't duplicate rows).
    combined_path = BASE / "docking_results_apo.csv"
    original = []
    with open(combined_path) as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames
        for r in reader:
            if r["mol_id"] not in LIGANDS:
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
