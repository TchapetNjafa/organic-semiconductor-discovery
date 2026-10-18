"""
Apo-receptor redocking of all 7 candidate molecules (JCIM desk-rejection fix, 2026-08-03).

Root cause fixed: `redock_candidates.py` docked the candidates into the ORIGINAL
receptor PDBQTs, which still contain the native co-crystallized ligand (holo),
while `calibrate_reference_ligands.py` stripped the native ligand to dock the
reference ligands into an apo receptor. The two protocols were therefore not
comparable -- exactly the objection raised by the JCIM associate editor
(reference ligands redocked against apo receptors, candidates docked against
holo receptors).

Fix: redock the same 7 candidates (1712, 17574, 18506, 9168 = stable donors;
4550, 17851, 20778 = reactive structures) against the SAME apo receptors
already generated for the reference-ligand calibration
(`calibration_outputs/{pdb_id}_apo.pdbqt`), using the identical protocol
(seed=42, exhaustiveness=32, num_modes=9) and the identical box (computed from
the full original PDB, unchanged by ligand stripping -- confirmed identical to
`calibration_reference_ligands.csv` box values).

This produces `docking_results_apo.csv`, directly comparable to
`calibration_reference_ligands.csv` (now both apo-vs-apo).
"""
import csv
from pathlib import Path

SMINA = "/home/tchapet/SMINA/smina.static"
BASE = Path(__file__).resolve().parent
LIGAND_DIR = BASE / "ligands"
APO_DIR = BASE / "calibration_outputs"
OUT_DIR = BASE / "outputs_apo"
LOG_DIR = BASE / "logs_apo"
OUT_DIR.mkdir(exist_ok=True)
LOG_DIR.mkdir(exist_ok=True)

SEED = 42
EXHAUSTIVENESS = 32
NUM_MODES = 9

MOL_IDS = ["1712", "17574", "18506", "9168", "4550", "17851", "20778"]

# target -> (pdb_id, box) -- box values taken verbatim from
# calibration_reference_ligands.csv (computed from the full original PDB,
# identical for holo and apo since stripping the ligand from the PDBQT does
# not change the full-structure heavy-atom extent used to define the box).
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


def run_smina(ligand_pdbqt: Path, receptor_pdbqt: Path, box: dict, out_pdbqt: Path, log_path: Path) -> float:
    import subprocess
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
    rows = []
    for target, (pdb_id, box) in TARGETS.items():
        apo_receptor = APO_DIR / f"{pdb_id}_apo.pdbqt"
        if not apo_receptor.exists():
            raise FileNotFoundError(f"missing apo receptor: {apo_receptor}")
        for mol_id in MOL_IDS:
            ligand_pdbqt = LIGAND_DIR / f"{mol_id}.pdbqt"
            if not ligand_pdbqt.exists():
                raise FileNotFoundError(f"missing ligand: {ligand_pdbqt}")
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

    results_path = BASE / "docking_results_apo.csv"
    with open(results_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {len(rows)} results to {results_path}")


if __name__ == "__main__":
    main()
