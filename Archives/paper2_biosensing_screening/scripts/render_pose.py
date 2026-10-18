"""
Story 1.4: render 1712's best docking pose (Hsp90/2XJX, -7.04 kcal/mol,
model 1 of docking_rerun/outputs/1712_Hsp90.pdbqt) as a publication figure.
Run with: pymol -cq render_pose.py
"""
from pymol import cmd

RECEPTOR = (
    "/home/tchapet/Post-Doc/ARTICLES-DOSSIERS-DES-ARTICLES-EN-REDACTION/"
    "NOUVELS-AXES-DE-RECHERCHE-A-REGARDER-URGEMMENT/ARTICLES-EN-REDACTIONS/"
    "ARTICLE_MVOTO/MVOTOOPV/PUBCHEMQC_MULTIFONCTIONAL_2025/receptors/2XJX.pdb"
)
POSE = (
    "/home/tchapet/Post-Doc/ARTICLES-DOSSIERS-DES-ARTICLES-EN-REDACTION/"
    "NOUVELS-AXES-DE-RECHERCHE-A-REGARDER-URGEMMENT/ARTICLES-EN-REDACTIONS/"
    "ARTICLE_MVOTO/Paper2_Sensors_ActuatorsB/docking_rerun/outputs/1712_Hsp90_pose.pdb"
)
OUT_PNG = (
    "/home/tchapet/Post-Doc/ARTICLES-DOSSIERS-DES-ARTICLES-EN-REDACTION/"
    "NOUVELS-AXES-DE-RECHERCHE-A-REGARDER-URGEMMENT/ARTICLES-EN-REDACTIONS/"
    "ARTICLE_MVOTO/Paper2_Sensors_ActuatorsB/figures/fig4_docking_pose.png"
)

cmd.load(RECEPTOR, "receptor")
cmd.load(POSE, "ligand")

cmd.bg_color("white")
cmd.hide("everything")

cmd.show("cartoon", "receptor")
cmd.color("gray80", "receptor")
cmd.set("cartoon_transparency", 0.2)

cmd.show("sticks", "ligand")
cmd.color("green", "ligand and elem C")
cmd.util.cnc("ligand")
cmd.set("stick_radius", 0.18, "ligand")

# Show only receptor residues near the ligand for a clean binding-site view.
cmd.select("pocket", "byres (receptor within 5 of ligand)")
cmd.show("sticks", "pocket and not name C+N+O")
cmd.color("gray60", "pocket and elem C")
cmd.util.cnc("pocket")

cmd.orient("ligand")
cmd.center("ligand or pocket")
cmd.zoom("ligand or pocket", 6)

cmd.set("ray_opaque_background", 1)
cmd.set("antialias", 2)
cmd.ray(2200, 2200)
cmd.png(OUT_PNG, dpi=600)
print(f"wrote {OUT_PNG}")
# Cropping to actual content happens in crop_pose.py (venv Python, has PIL/numpy
# -- PyMOL's own bundled Python does not).
