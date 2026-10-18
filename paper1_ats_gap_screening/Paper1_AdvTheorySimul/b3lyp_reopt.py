"""B3LYP/6-31G* geometry optimisation (geomeTRIC) then KS gap; closes the geometry caveat of L-2026-10-17-03.

usage: python -I b3lyp_reopt.py <in.xyz> <out_dir>
"""
import os
import sys

from pyscf import dft, gto, lib
from pyscf.geomopt.geometric_solver import kernel as geomopt

HARTREE_EV = 27.211386245988
lib.num_threads(8)
src, out = sys.argv[1], sys.argv[2]
tag = f"{os.path.basename(os.path.dirname(src))}_{os.path.splitext(os.path.basename(src))[0]}"
mf = dft.RKS(gto.M(atom=src, basis="6-31g*", verbose=0), xc="b3lyp")
conv_geom, mol_eq = geomopt(mf, maxsteps=150)
if not conv_geom:
    sys.exit(f"{tag}: geometry optimisation did not converge")
mol_eq.tofile(os.path.join(out, f"{tag}_b3lyp.xyz"))
rks = dft.RKS(mol_eq, xc="b3lyp").run()
occ, e = rks.mo_occ, rks.mo_energy
homo, lumo = e[occ > 0].max() * HARTREE_EV, e[occ == 0].min() * HARTREE_EV
print(f"{tag}\tgeom_converged={conv_geom}\tconverged={rks.converged}\tHOMO={homo:.3f}\tLUMO={lumo:.3f}\tGAP={lumo - homo:.3f} eV", flush=True)
