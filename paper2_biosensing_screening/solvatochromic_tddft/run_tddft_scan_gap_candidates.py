"""
TD-DFT solvatochromic scan for the 3 heavy-atom-gap-filling candidates
(20549, 5081, 23424; JCIM desk-rejection fix #4, 2026-08-03).

Identical protocol and rationale to run_tddft_scan.py, applied to the same
optimized geometries (geometries/{mol_id}.xyz, B3LYP/6-31G*-optimized by
optimize_gap_geometries.py) and the same functional pair (B3LYP -- matching
Paper 1's DFT level -- and CAM-B3LYP -- range-separated cross-check for
charge-transfer character), same two solvents (Toluene = nonpolar OPV-blend
proxy, Water = polar biosensing-interface proxy), same TDA/8-root/CPCM
settings. Kept as a separate script (rather than extending
run_tddft_scan.py's MOLECULES list) so re-running it does not also re-run
the already-completed TD-DFT for the original 7 candidates.

Usage:
    python3 run_tddft_scan_gap_candidates.py                       # B3LYP, run
    python3 run_tddft_scan_gap_candidates.py --functional CAM-B3LYP
    python3 run_tddft_scan_gap_candidates.py --parse-only           # re-parse only
"""
import argparse
import os
import re
import subprocess
from pathlib import Path

HOME = Path.home()
ORCA_EXE = HOME / "orca-6-1-1" / "orca"

BASE = Path(__file__).resolve().parent
GEOM_DIR = BASE / "geometries"

MOLECULES = ["20549", "5081", "23424"]
CHARGE = 0
MULTIPLICITY = 1

SOLVENTS = ["Toluene", "Water"]
NROOTS = 8


def _physical_cores() -> int:
    """Open MPI counts physical cores as slots; logical (SMT) count over-
    subscribes and aborts at startup. Derive physical cores from core ids."""
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

FUNCTIONAL_KEYWORDS = {
    "B3LYP": "B3LYP 6-31G*",
    "CAM-B3LYP": "CAM-B3LYP 6-31G*",
}


def write_input(mol_id: str, solvent: str, nprocs: int, functional: str, run_dir: Path) -> Path:
    run_subdir = run_dir / f"{mol_id}_{solvent}"
    run_subdir.mkdir(parents=True, exist_ok=True)
    geom_src = GEOM_DIR / f"{mol_id}.xyz"
    geom_dst = run_subdir / f"{mol_id}.xyz"
    geom_dst.write_text(geom_src.read_text())

    inp_path = run_subdir / f"{mol_id}_{solvent}.inp"
    inp_path.write_text(
        f"! {FUNCTIONAL_KEYWORDS[functional]} CPCM({solvent}) TightSCF\n"
        f"%pal nprocs {nprocs} end\n"
        f"%maxcore {MAXCORE_MB}\n"
        f"%tddft\n"
        f"  NRoots {NROOTS}\n"
        f"  TDA true\n"
        f"end\n"
        f"* xyzfile {CHARGE} {MULTIPLICITY} {geom_dst.name}\n"
    )
    return inp_path


def run_orca(inp_path: Path) -> Path:
    out_path = inp_path.with_suffix(".out")
    with open(out_path, "w") as fh:
        subprocess.run(
            [str(ORCA_EXE), str(inp_path)],
            stdout=fh,
            stderr=subprocess.STDOUT,
            cwd=inp_path.parent,
            check=False,
        )
    return out_path


ABS_TABLE_ROW = re.compile(
    r"^\s*\d+-\d[A-Za-z\"']*\s*->\s*(\d+)-\d[A-Za-z\"']*\s+"
    r"([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)"
)


def parse_states(out_path: Path):
    text = out_path.read_text(errors="replace")
    if "ORCA TERMINATED NORMALLY" not in text:
        return None

    lines = text.splitlines()
    states = []
    in_table = False
    for line in lines:
        if "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in line:
            in_table = True
            continue
        if in_table:
            m = ABS_TABLE_ROW.match(line)
            if m:
                state, energy_eV, _energy_cm1, wavelength_nm, fosc = m.groups()
                states.append({
                    "state": int(state),
                    "energy_eV": float(energy_eV),
                    "wavelength_nm": float(wavelength_nm),
                    "oscillator_strength": float(fosc),
                })
            elif states and not line.strip():
                break
    return states or None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--parse-only", action="store_true")
    parser.add_argument("--nprocs", type=int, default=DEFAULT_NPROCS)
    parser.add_argument("--functional", choices=list(FUNCTIONAL_KEYWORDS), default="B3LYP")
    args = parser.parse_args()

    run_dir = BASE / f"runs_gap_{args.functional}"
    results_csv = BASE / f"solvatochromic_results_gap_candidates_{args.functional}.csv"
    run_dir.mkdir(exist_ok=True)
    rows = []

    for mol_id in MOLECULES:
        for solvent in SOLVENTS:
            inp_path = write_input(mol_id, solvent, args.nprocs, args.functional, run_dir)
            out_path = inp_path.with_suffix(".out")

            if not args.parse_only:
                if not ORCA_EXE.exists():
                    raise FileNotFoundError(f"ORCA not found at {ORCA_EXE}")
                print(f"Running {mol_id} / {solvent} ({args.functional}) ...")
                run_orca(inp_path)

            states = parse_states(out_path)
            if states is None:
                print(f"  WARNING: {mol_id}/{solvent} did not complete normally "
                      f"or produced no absorption table -- see {out_path}")
                continue
            for s in states:
                rows.append({"mol_id": mol_id, "solvent": solvent, **s})
            print(f"  {mol_id}/{solvent}: S1 = {states[0]['energy_eV']:.3f} eV "
                  f"({states[0]['wavelength_nm']:.1f} nm), "
                  f"f = {states[0]['oscillator_strength']:.4f}")

    if not rows:
        print("No results parsed -- nothing written.")
        return

    import csv
    with open(results_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {len(rows)} state rows to {results_csv}")

    print("\nSolvatochromic shift (S1, water - toluene):")
    s1 = {(r["mol_id"], r["solvent"]): r for r in rows if r["state"] == 1}
    for mol_id in MOLECULES:
        tol = s1.get((mol_id, "Toluene"))
        wat = s1.get((mol_id, "Water"))
        if tol and wat:
            shift_eV = wat["energy_eV"] - tol["energy_eV"]
            shift_nm = wat["wavelength_nm"] - tol["wavelength_nm"]
            print(f"  {mol_id}: {shift_eV:+.4f} eV ({shift_nm:+.2f} nm)")


if __name__ == "__main__":
    main()
