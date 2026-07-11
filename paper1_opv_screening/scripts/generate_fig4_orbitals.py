#!/usr/bin/env python3
"""
Figure 4: real HOMO/LUMO isosurfaces for the paired case studies
17851 (top-ranked, unstable azido-triazine) and 1712 (the stable survivor),
rendered directly from the pyscf .molden wavefunctions in NTO_and_related_study.
No new QM: orbitals come from existing files.

2 rows (molecule) x 2 cols (HOMO, LUMO). Positive lobe red, negative lobe blue.
Emits figures/figure4_orbitals.{pdf,png}.
"""
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
import figstyle as fs
from skimage import measure
from pyscf.tools import molden


MOL = {
    17851: "17851  (azido-triazine, unstable)",
    1712: "1712  (stable acceptor)",
}
BOHR = 0.529177
ELEM_COLOR = {1: "#e0e0e0", 6: "#404040", 7: "#2b6cb0", 8: "#c53030", 16: "#d69e2e"}
ELEM_SIZE = {1: 10, 6: 26, 7: 30, 8: 30, 16: 44}


def load(mid):
    res = molden.load(f"../../data/data/nbo_analysis/{mid}/{mid}.molden")
    mol, mo_e, mo_c, mo_occ = res[0], res[1], res[2], np.asarray(res[3])
    homo = np.where(mo_occ > 0)[0][-1]
    return mol, mo_c, homo, homo + 1


def orbital_volume(mol, coeff, n=76, pad=3.2):
    xyz = mol.atom_coords()  # Bohr
    lo, hi = xyz.min(0) - pad, xyz.max(0) + pad
    axes = [np.linspace(lo[i], hi[i], n) for i in range(3)]
    grid = np.array(np.meshgrid(*axes, indexing="ij")).reshape(3, -1).T
    ao = mol.eval_gto("GTOval", grid)
    vol = (ao @ coeff).reshape(n, n, n)
    spacing = tuple((hi - lo) / (n - 1))
    return vol, lo, spacing, xyz


def draw_orbital(ax, mol, coeff, title):
    vol, lo, spacing, xyz = orbital_volume(mol, coeff)
    iso = 0.035
    m = np.abs(vol).max()
    if iso > 0.9 * m:                     # fall back if orbital is diffuse
        iso = 0.35 * m
    for sign, color in ((+1, "#e05252"), (-1, "#5b8dd6")):
        data = sign * vol
        if data.max() <= iso:
            continue
        verts, faces, _, _ = measure.marching_cubes(data, level=iso, spacing=spacing)
        verts = (verts + lo) * BOHR       # Bohr -> Angstrom
        ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2],
                        color=color, alpha=0.55, linewidth=0, antialiased=True, shade=True)
    a = xyz * BOHR
    for i in range(mol.natm):
        z = mol.atom_charge(i)
        ax.scatter(*a[i], s=ELEM_SIZE.get(z, 20), c=ELEM_COLOR.get(z, "#7f7f7f"),
                   edgecolors="black", linewidths=0.3, depthshade=True, zorder=5)
    # bonds
    for i in range(mol.natm):
        for j in range(i + 1, mol.natm):
            d = np.linalg.norm(a[i] - a[j])
            if d < 1.75:
                ax.plot(*zip(a[i], a[j]), color="#2d3748", linewidth=1.0, zorder=4)
    ax.set_title(title, fontsize=8.5, pad=-2)
    ax.set_box_aspect((1, 1, 1), zoom=1.7)
    ax.set_axis_off()
    c = a.mean(0); r = (a.max(0) - a.min(0)).max() / 2 + 0.8
    ax.set_xlim(c[0]-r, c[0]+r); ax.set_ylim(c[1]-r, c[1]+r); ax.set_zlim(c[2]-r, c[2]+r)
    ax.view_init(elev=22, azim=-60)


fig = plt.figure(figsize=(fs.DOUBLE, fs.DOUBLE*0.68))
for row, mid in enumerate(MOL):
    mol, mo_c, homo, lumo = load(mid)
    for col, (idx, name) in enumerate([(homo, "HOMO"), (lumo, "LUMO")]):
        ax = fig.add_subplot(2, 2, row * 2 + col + 1, projection="3d")
        draw_orbital(ax, mol, mo_c[:, idx], f"Molecule {mid} — {name}")
fig.text(0.5, 0.975, "HOMO and LUMO isosurfaces (isovalue 0.035 a.u.; red $+$, blue $-$)",
         ha="center", fontsize=9)
fig.subplots_adjust(left=0.0, right=1.0, top=0.94, bottom=0.0, wspace=0.0, hspace=0.02)
for ext in ("pdf", "png"):
    fig.savefig(f"figures/figure4_orbitals.{ext}", dpi=500, bbox_inches="tight")
print("wrote figures/figure4_orbitals.pdf / .png")
