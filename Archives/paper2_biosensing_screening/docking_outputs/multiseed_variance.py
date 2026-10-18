"""
Multi-seed docking variance study (addresses Phase-5 review item H5).

The original screen (docking_results_clean.csv) used a single smina seed (42).
A reviewer cannot judge whether the reported affinity differences exceed the
docking search's own run-to-run scatter. This script re-docks every
ligand x target pair under THREE seeds, holding everything else fixed
(same prepared ligand PDBQT conformer, same blind-docking box, same
exhaustiveness/num_modes), and reports per-pair mean +/- SD of the best-mode
affinity.

Isolating only --seed is deliberate: it measures the stochastic-search
component of the score, which is the quantity the single-seed protocol left
unquantified. Ligand-conformer and box definitions are held identical to the
published run so the seed-42 column reproduces docking_results_clean.csv.

Idempotent: per (mol, target, seed) results are cached under logs_multiseed/;
re-running skips completed cells and only fills gaps.
"""
import csv
import statistics
import subprocess
from pathlib import Path

SMINA = "/home/tchapet/SMINA/smina.static"
BASE = Path(__file__).resolve().parent
RECEPTOR_DIR = Path(
    "/home/tchapet/Post-Doc/ARTICLES-DOSSIERS-DES-ARTICLES-EN-REDACTION/"
    "NOUVELS-AXES-DE-RECHERCHE-A-REGARDER-URGEMMENT/ARTICLES-EN-REDACTIONS/"
    "ARTICLE_MVOTO/MVOTOOPV/PUBCHEMQC_MULTIFONCTIONAL_2025/receptors"
)
LIGAND_DIR = BASE / "ligands"
OUT_DIR = BASE / "outputs_multiseed"
LOG_DIR = BASE / "logs_multiseed"
OUT_DIR.mkdir(exist_ok=True)
LOG_DIR.mkdir(exist_ok=True)

SEEDS = [42, 7, 123]
EXHAUSTIVENESS = 32
NUM_MODES = 9

RECEPTORS = {
    "HIV1protease": "1DMP",
    "Hsp90": "2XJX",
    "Neurodegenerative_1SYH": "1SYH",
    "COVID19_6Y2F": "6Y2F",
}


def load_boxes(clean_csv: Path) -> dict:
    """Read the exact per-target box used in the published single-seed run."""
    boxes = {}
    with open(clean_csv) as fh:
        for r in csv.DictReader(fh):
            t = r["target"]
            if t not in boxes:
                boxes[t] = {
                    "center": (float(r["box_center_x"]), float(r["box_center_y"]), float(r["box_center_z"])),
                    "size": (float(r["box_size_x"]), float(r["box_size_y"]), float(r["box_size_z"])),
                }
    return boxes


def parse_best_affinity(log_path: Path) -> float:
    with open(log_path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 2 and parts[0] == "1":
                try:
                    return float(parts[1])
                except ValueError:
                    continue
    raise RuntimeError(f"could not parse affinity from {log_path}")


def run_one(mol_id: str, target: str, pdb_id: str, seed: int, box: dict) -> float:
    log_path = LOG_DIR / f"{mol_id}_{target}_seed{seed}.log"
    out_pdbqt = OUT_DIR / f"{mol_id}_{target}_seed{seed}.pdbqt"
    if log_path.exists():
        try:
            return parse_best_affinity(log_path)
        except RuntimeError:
            pass  # incomplete/corrupt log -> re-run
    cx, cy, cz = box["center"]
    sx, sy, sz = box["size"]
    cmd = [
        SMINA,
        "--receptor", str(RECEPTOR_DIR / f"{pdb_id}.pdbqt"),
        "--ligand", str(LIGAND_DIR / f"{mol_id}.pdbqt"),
        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
        "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", str(NUM_MODES),
        "--seed", str(seed),
        "--out", str(out_pdbqt),
        "--log", str(log_path),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"smina failed {mol_id}/{target}/seed{seed}: {result.stderr}")
    return parse_best_affinity(log_path)


def main():
    boxes = load_boxes(BASE / "docking_results_clean.csv")
    ligands = sorted(p.stem for p in LIGAND_DIR.glob("*.pdbqt"))
    rows = []
    for target, pdb_id in RECEPTORS.items():
        box = boxes[target]
        for mol_id in ligands:
            affs = []
            for seed in SEEDS:
                a = run_one(mol_id, target, pdb_id, seed, box)
                affs.append(a)
                print(f"{mol_id} {target} seed{seed}: {a}", flush=True)
            mean = statistics.mean(affs)
            sd = statistics.stdev(affs) if len(affs) > 1 else 0.0
            row = {"mol_id": mol_id, "target": target, "pdb_id": pdb_id,
                   "exhaustiveness": EXHAUSTIVENESS, "num_modes": NUM_MODES}
            for seed, a in zip(SEEDS, affs):
                row[f"aff_seed{seed}"] = a
            row["mean_kcal_mol"] = round(mean, 3)
            row["sd_kcal_mol"] = round(sd, 3)
            row["range_kcal_mol"] = round(max(affs) - min(affs), 3)
            rows.append(row)

    out_csv = BASE / "docking_multiseed_variance.csv"
    with open(out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    max_sd = max(r["sd_kcal_mol"] for r in rows)
    mean_sd = statistics.mean(r["sd_kcal_mol"] for r in rows)
    print(f"\nWrote {len(rows)} pairs to {out_csv}")
    print(f"SD across all pairs: mean={mean_sd:.3f}, max={max_sd:.3f} kcal/mol")


if __name__ == "__main__":
    main()
