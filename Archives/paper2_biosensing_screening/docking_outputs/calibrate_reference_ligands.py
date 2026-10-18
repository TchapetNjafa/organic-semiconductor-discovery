"""
Reference-ligand calibration for the four docking targets (JCIM revision, 2026-07-14).

Purpose: answer the "uncalibrated docking" reviewer objection by redocking each
target's native, co-crystallized ligand under the IDENTICAL protocol used for the
four organic-semiconductor candidates (same whole-protein box, seed=42,
exhaustiveness=32, num_modes=9). The best-mode affinity of a known binder places
the candidate affinities on an interpretable scale.

Native ligands (resname in the deposited crystal structures):
  1DMP (HIV-1 protease)      : DMQ  (cyclic-urea inhibitor, XK263/DMP323 family)
  2XJX (Hsp90)               : XJX  (co-crystallized inhibitor)
  1SYH (ionotropic GluR2)    : CPW  (co-crystallized antagonist)
  6Y2F (SARS-CoV-2 M-pro)    : O6K  (alpha-ketoamide 13b)

Protocol notes / honesty:
- The receptor PDBQT files used for the candidate runs still CONTAIN their native
  ligand. For a valid redock the native ligand is stripped from the receptor
  (-> apo receptor) here; the box is recomputed identically from the original full
  PDB so it matches the candidate runs bit-for-bit.
- The native crystal ligand is docked from its crystallographic coordinates with
  hydrogens added by Open Babel and Gasteiger charges assigned (crystal ligands
  lack H). No RDKit re-embedding is applied (the crystal conformer is the input),
  which is standard redocking practice and is stated as such.
- We report the reference best-mode affinity per target and, where obrms is
  available, the RMSD of the best redocked pose to the crystal pose (a pose-recovery
  check). RMSD < ~2 A is the usual redocking-success threshold.
"""
import csv
import shutil
import subprocess
from pathlib import Path

SMINA = "/home/tchapet/SMINA/smina.static"
BASE = Path(__file__).resolve().parent
RECEPTOR_DIR = Path(
    "/home/tchapet/Post-Doc/ARTICLES-DOSSIERS-DES-ARTICLES-EN-REDACTION/"
    "NOUVELS-AXES-DE-RECHERCHE-A-REGARDER-URGEMMENT/ARTICLES-EN-REDACTIONS/"
    "ARTICLE_MVOTO/MVOTOOPV/PUBCHEMQC_MULTIFONCTIONAL_2025/receptors"
)
OUT_DIR = BASE / "calibration_outputs"
OUT_DIR.mkdir(exist_ok=True)

SEED = 42
EXHAUSTIVENESS = 32
NUM_MODES = 9
PADDING_ANGSTROM = 8.0

# target -> (pdb_id, native-ligand resname)
TARGETS = {
    "HIV1protease": ("1DMP", "DMQ"),
    "Hsp90": ("2XJX", "XJX"),
    "Neurodegenerative_1SYH": ("1SYH", "CPW"),
    "COVID19_6Y2F": ("6Y2F", "O6K"),
}
WATERS = {"HOH", "WAT", "DOD"}


def receptor_box(receptor_pdb: Path) -> dict:
    """Identical box computation to redock_candidates.py (full-structure extent)."""
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
    center = {"x": (max(xs) + min(xs)) / 2,
              "y": (max(ys) + min(ys)) / 2,
              "z": (max(zs) + min(zs)) / 2}
    size = {"x": (max(xs) - min(xs)) + 2 * PADDING_ANGSTROM,
            "y": (max(ys) - min(ys)) + 2 * PADDING_ANGSTROM,
            "z": (max(zs) - min(zs)) + 2 * PADDING_ANGSTROM}
    return {"center": center, "size": size}


def extract_native_ligand(pdb: Path, resname: str, out_pdb: Path) -> int:
    """Write only the HETATM records of `resname` (crystal ligand) to out_pdb."""
    n = 0
    with open(pdb) as fh, open(out_pdb, "w") as out:
        for line in fh:
            if line.startswith("HETATM") and line[17:20].strip() == resname:
                out.write(line)
                n += 1
        out.write("END\n")
    return n


def strip_ligand_from_receptor(receptor_pdbqt: Path, resname: str, out_pdbqt: Path) -> int:
    """Copy receptor PDBQT minus the coordinate records of `resname` (-> apo)."""
    removed = 0
    with open(receptor_pdbqt) as fh, open(out_pdbqt, "w") as out:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() == resname:
                removed += 1
                continue
            out.write(line)
    return removed


def prep_ligand_pdbqt(lig_pdb: Path, lig_pdbqt: Path) -> None:
    subprocess.run(
        ["obabel", str(lig_pdb), "-O", str(lig_pdbqt),
         "-h", "--partialcharge", "gasteiger"],
        check=True, capture_output=True, text=True,
    )


def run_smina(ligand_pdbqt: Path, receptor_pdbqt: Path, box: dict,
              out_pdbqt: Path, log_path: Path) -> float:
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
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"smina failed: {res.stderr}\n{res.stdout}")
    with open(log_path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 2 and parts[0] == "1":
                try:
                    return float(parts[1])
                except ValueError:
                    continue
    raise RuntimeError(f"could not parse best-mode affinity from {log_path}")


def try_rmsd(crystal_pdb: Path, docked_pdbqt: Path) -> str:
    """Best-effort pose-recovery RMSD via obrms; returns '' if unavailable."""
    if shutil.which("obrms") is None:
        return ""
    # Convert docked best pose (model 1) to pdb for obrms.
    best_pdb = docked_pdbqt.with_suffix(".best.pdb")
    try:
        subprocess.run(["obabel", str(docked_pdbqt), "-O", str(best_pdb), "-l", "1"],
                       check=True, capture_output=True, text=True)
        res = subprocess.run(["obrms", str(crystal_pdb), str(best_pdb)],
                             capture_output=True, text=True)
        out = res.stdout.strip()
        return out.split()[-1] if out else ""
    except Exception:
        return ""


def main():
    rows = []
    for target, (pdb_id, resname) in TARGETS.items():
        full_pdb = RECEPTOR_DIR / f"{pdb_id}.pdb"
        rec_pdbqt = RECEPTOR_DIR / f"{pdb_id}.pdbqt"
        box = receptor_box(full_pdb)

        lig_pdb = OUT_DIR / f"{pdb_id}_{resname}_native.pdb"
        lig_pdbqt = OUT_DIR / f"{pdb_id}_{resname}_native.pdbqt"
        apo_pdbqt = OUT_DIR / f"{pdb_id}_apo.pdbqt"
        out_pose = OUT_DIR / f"{pdb_id}_{resname}_redock.pdbqt"
        log_path = OUT_DIR / f"{pdb_id}_{resname}_redock.log"

        n_lig = extract_native_ligand(full_pdb, resname, lig_pdb)
        n_removed = strip_ligand_from_receptor(rec_pdbqt, resname, apo_pdbqt)
        prep_ligand_pdbqt(lig_pdb, lig_pdbqt)
        aff = run_smina(lig_pdbqt, apo_pdbqt, box, out_pose, log_path)
        rmsd = try_rmsd(lig_pdb, out_pose)

        print(f"{pdb_id} {resname}: ref affinity = {aff} kcal/mol "
              f"(ligand atoms={n_lig}, removed from receptor={n_removed}, RMSD={rmsd or 'n/a'})")
        rows.append({
            "target": target, "pdb_id": pdb_id, "native_ligand": resname,
            "ref_affinity_kcal_mol": aff, "redock_rmsd_angstrom": rmsd,
            "n_ligand_atoms": n_lig, "n_removed_from_receptor": n_removed,
            "box_center_x": box["center"]["x"], "box_center_y": box["center"]["y"],
            "box_center_z": box["center"]["z"], "box_size_x": box["size"]["x"],
            "box_size_y": box["size"]["y"], "box_size_z": box["size"]["z"],
            "exhaustiveness": EXHAUSTIVENESS, "seed": SEED, "num_modes": NUM_MODES,
        })

    out_csv = BASE / "calibration_reference_ligands.csv"
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {len(rows)} reference-ligand results to {out_csv}")


if __name__ == "__main__":
    main()
