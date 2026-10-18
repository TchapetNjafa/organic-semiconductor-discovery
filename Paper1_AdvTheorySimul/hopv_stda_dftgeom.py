"""Same sTDA-xTB step as hopv_stda.py, but on HOPV15's own DFT conformer-1 geometry (isolates the UFF-geometry effect).

usage: python -I hopv_stda_dftgeom.py <HOPV_15_revised_2.data> <out_dir> [n_jobs]
writes <out_dir>/hopv15_stda_dftgeom.csv (small molecules only: construction=molecule, complement != TiO2)
"""
import csv
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from hopv_stda import EXP_COLS, STDA, XTB4STDA


def conformer1(path):
    lines = Path(path).read_text().splitlines()
    i, k = 0, 0
    while i < len(lines):
        exp = dict(zip(EXP_COLS, (x.strip() for x in lines[i + 2].split(","))))
        n_conf = int(lines[i + 4])
        i += 5
        for c in range(n_conf):
            n_at = int(lines[i + 1])
            if c == 0:
                xyz = f"{n_at}\nhopv {k}\n" + "\n".join(lines[i + 2: i + 2 + n_at]) + "\n"
            i += 2 + n_at + 4
        yield k, exp, xyz
        k += 1


def run(job):
    idx, xyz, out_dir = job
    wd = Path(out_dir) / f"dft_{idx:03d}"
    wd.mkdir(parents=True, exist_ok=True)
    (wd / "mol.xyz").write_text(xyz)
    env = {**os.environ, "XTB4STDAHOME": os.path.expanduser("~/xtb4stda"), "OMP_NUM_THREADS": "1"}
    for stale in ("wfn.xtb", "tda.dat"):
        (wd / stale).unlink(missing_ok=True)
    a = subprocess.run([XTB4STDA, "mol.xyz"], cwd=wd, env=env, capture_output=True, timeout=1800)
    b = subprocess.run([STDA, "-xtb", "wfn.xtb"], cwd=wd, env=env, capture_output=True, timeout=1800)
    if a.returncode or b.returncode:
        return (idx, "nonzero_exit", None, None)
    tda = wd / "tda.dat"
    if not tda.exists():
        return idx, "stda_failed", None, None
    (wd / "wfn.xtb").unlink(missing_ok=True)
    rows = [l.split() for l in tda.read_text().split("DATXY", 1)[1].splitlines() if len(l.split()) >= 3]
    return (idx, "success", float(rows[0][1]), float(rows[0][2])) if rows else (idx, "no_states", None, None)


if __name__ == "__main__":
    src, out = sys.argv[1], sys.argv[2]
    jobs = [(k, xyz, out) for k, e, xyz in conformer1(src)
            if e["construction"] == "molecule" and "TiO2" not in e["complement"]]
    assert len(jobs) == 63, len(jobs)
    with ProcessPoolExecutor(int(sys.argv[3]) if len(sys.argv) > 3 else 6) as ex, \
            open(Path(out) / "hopv15_stda_dftgeom.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["idx", "status", "E_S1_dftgeom_eV", "f_S1_dftgeom"])
        for row in ex.map(run, jobs):
            w.writerow(row)
            f.flush()
