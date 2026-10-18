"""B3LYP/6-31G* (Cartesian) ground state with PySCF -> Molden -> std2 sTDA and XsTDA lowest singlet.

usage: python -I std2_b3lyp.py <out_dir> <xyz> [<xyz> ...]   (env NTHREADS, default 8)
prints: tag  KSgap  S1_sTDA f  S1_XsTDA f   (eV)
Notes: std2 needs "[Atoms] Angs" (it mis-reads PySCF's "(AU)" header) and raw PySCF Cartesian coefficients.
std2 binary needs MKL; LD_LIBRARY_PATH must contain libmkl_intel_ilp64.so.
"""
import os
import re
import subprocess
import sys
from pathlib import Path

from pyscf import dft, gto, lib
from pyscf.tools import molden

BOHR = 0.52917721092
HARTREE_EV = 27.211386245988
STD2 = os.path.expanduser("~/std2-2.0.1/_build/std2")
lib.num_threads(int(os.environ.get("NTHREADS", "8")))
STATE = re.compile(r"^\s+1\s+([0-9.]+)\s+[0-9.]+\s+([0-9.]+)\s")


def to_angs_header(src, dst):
    out, in_atoms = [], False
    for line in Path(src).read_text().splitlines():
        if line.startswith("[Atoms]"):
            out.append("[Atoms] Angs")
            in_atoms = True
            continue
        if line.startswith("["):
            in_atoms = False
        f = line.split()
        if in_atoms and len(f) == 6:
            line = f"{f[0]} {f[1]} {f[2]} " + " ".join(f"{float(x) * BOHR:.10f}" for x in f[3:])
        out.append(line)
    Path(dst).write_text("\n".join(out) + "\n")


def s1(out_text):
    """First line of std2's final state table: energy (eV) and f."""
    if "excitation energies, transition moments" not in out_text:
        return float("nan"), float("nan")
    block = out_text.split("excitation energies, transition moments", 1)[1]
    for line in block.splitlines():
        m = STATE.match(line)
        if m:
            return float(m[1]), float(m[2])
    return float("nan"), float("nan")


def run(xyz, out_dir):
    tag = f"{Path(xyz).parent.name}_{Path(xyz).stem}"
    wd = Path(out_dir) / tag
    wd.mkdir(parents=True, exist_ok=True)
    mol = gto.M(atom=xyz, basis="6-31g*", cart=True, verbose=0)
    mf = dft.RKS(mol, xc="b3lyp").run()
    if not mf.converged:
        print(f"{tag}\tconv=False\tSKIPPED (unconverged SCF)", flush=True)
        return
    occ, e = mf.mo_occ, mf.mo_energy
    gap = (e[occ == 0].min() - e[occ > 0].max()) * HARTREE_EV
    molden.from_scf(mf, str(wd / "raw.molden"))
    to_angs_header(wd / "raw.molden", wd / "molden.input")
    env = {**os.environ, "OMP_NUM_THREADS": os.environ.get("NTHREADS", "8"),
           "MKL_NUM_THREADS": os.environ.get("NTHREADS", "8")}
    res = []
    for extra, name in (([], "stda"), (["-XsTD"], "xstda")):
        p = subprocess.run([STD2, "-f", "molden.input", "-sty", "1", "-ax", "0.2", "-e", "7", *extra],
                           cwd=wd, env=env, capture_output=True, text=True, timeout=7200)
        (wd / f"{name}.out").write_text(p.stdout + p.stderr)
        res.extend(s1(p.stdout) if p.returncode == 0 else (float("nan"), float("nan")))
    print(f"{tag}\tconv={mf.converged}\tKS={gap:.3f}\tS1_sTDA={res[0]:.3f}\tf={res[1]:.4f}"
          f"\tS1_XsTDA={res[2]:.3f}\tf={res[3]:.4f}", flush=True)


if __name__ == "__main__":
    for x in sys.argv[2:]:
        run(x, sys.argv[1])
