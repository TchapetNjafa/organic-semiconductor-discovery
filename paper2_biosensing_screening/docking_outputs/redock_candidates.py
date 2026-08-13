"""
Clean, reproducible re-docking of the 4 genuine Paper-1 OPV candidates
(1712, 17851, 20778, 4550) against 4 tractable single-domain receptors
(HIV1protease/1DMP, Hsp90/2XJX, Neurodegenerative/1SYH, COVID19/6Y2F).

Fixes applied vs the prior (unreliable) runs:
- smina called by full path (was silently failing when unqualified).
- Single deterministic 3D embedding (RDKit ETKDGv3, fixed seed) -> PDB ->
  PDBQT via obabel with --partialcharge gasteiger and NO further 3D
  generation (the prior pipeline re-ran OpenBabel's own unseeded make3D
  on top of the RDKit conformer, which is the likely source of
  run-to-run drift).
- Whole-protein blind-docking box computed per receptor from heavy-atom
  extent + fixed padding (documented, reproducible), not a hardcoded
  (0,0,0)/20x20x20 box.
- Proteasome (7PG9, 28-chain/45k-atom 20S core particle) and the two
  broken receptors (4LDE, 5JWA) are excluded -- not in scope here.
- Fixed exhaustiveness, fixed seed, fixed num_modes for every run.
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

LIGANDS = {
    "1712": "O=C1N=c2c(=C1CNc1ccc(cc1)S(=O)(=O)Nc1ccccn1)c1scnc1cc2",
    "4550": "O=NC1=c2cc(N(O)O)c3c(c2=NC1=O)CCCC3",
    "17851": "CC[N-]c1nc(N[N+]#N)nc(n1)NC(C)C",
    "20778": "N#[N+]Nc1nc([N-]C(C)C)nc(n1)SC",
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

    # Single conversion pass: PDB (already 3D) -> PDBQT. No --gen3d, so
    # OpenBabel does not regenerate/overwrite the RDKit conformer.
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

    results_path = BASE / "docking_results_clean.csv"
    with open(results_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {len(rows)} results to {results_path}")


if __name__ == "__main__":
    main()
