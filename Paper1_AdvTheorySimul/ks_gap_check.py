"""B3LYP/6-31G* KS gap on the xTB geometry, to test whether a PubChemQC gap is reproducible.

usage: python -I ks_gap_check.py <xyz> [<xyz> ...]
Closed-shell RKS plus a UKS triplet probe; prints one line per molecule.
"""
import os
import sys

from pyscf import dft, gto, lib

HARTREE_EV = 27.211386245988
lib.num_threads(int(os.environ.get("NTHREADS", "8")))


def gap(mf):
    occ, e = mf.mo_occ, mf.mo_energy
    homo, lumo = e[occ > 0].max(), e[occ == 0].min()
    return homo * HARTREE_EV, lumo * HARTREE_EV, (lumo - homo) * HARTREE_EV


for path in sys.argv[1:]:
    mol = gto.M(atom=path, basis=os.environ.get("BASIS", "6-31g*"), verbose=0)
    rks = dft.RKS(mol, xc="b3lyp").run()
    h, l, g = gap(rks)
    dE = float("nan")
    if os.environ.get("TRIPLET", "1") == "1":
        trip = gto.M(atom=path, basis=os.environ.get("BASIS", "6-31g*"), spin=2, verbose=0)
        uks = dft.UKS(trip, xc="b3lyp").run()
        dE = (uks.e_tot - rks.e_tot) * HARTREE_EV if (uks.converged and rks.converged) else float("nan")
    print(f"{path}\tconverged={rks.converged}\tHOMO={h:.3f}\tLUMO={l:.3f}\tGAP={g:.3f} eV\tE(T1)-E(S0)={dE:.3f} eV", flush=True)
