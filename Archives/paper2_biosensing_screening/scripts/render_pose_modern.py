"""
Fig 4 (modern restyle) -- 1712's best docking pose against Hsp90 (2XJX,
-7.0 kcal/mol, model 1). Restyled 2026-07-17 to match the modern-minimal
figure system: ligand in the stable-donor BLUE (1712 is a stable donor),
clean white cartoon, wider framing so the ligand is not clipped.
Run with: pymol -cq render_pose_modern.py
"""
from pymol import cmd

BASE = ("/home/tchapet/Post-Doc/ARTICLES-DOSSIERS-DES-ARTICLES-EN-REDACTION/"
        "NOUVELS-AXES-DE-RECHERCHE-A-REGARDER-URGEMMENT/ARTICLES-EN-REDACTIONS/"
        "ARTICLE_MVOTO/")
RECEPTOR = BASE + "MVOTOOPV/PUBCHEMQC_MULTIFONCTIONAL_2025/receptors/2XJX.pdb"
POSE = BASE + "Paper2_Sensors_ActuatorsB/docking_rerun/outputs/1712_Hsp90_pose.pdb"
OUT_PNG = BASE + "Paper2_JCIM/figures/fig4_docking_pose.png"

cmd.load(RECEPTOR, "receptor")
cmd.load(POSE, "ligand")

cmd.bg_color("white")
cmd.hide("everything")

# Receptor: clean light cartoon, faint so the ligand reads first.
cmd.show("cartoon", "receptor")
cmd.color("gray90", "receptor")
cmd.set("cartoon_transparency", 0.35, "receptor")
cmd.set("cartoon_fancy_helices", 1)

# Pocket surface hint: soft translucent surface on nearby residues only.
cmd.select("pocket", "byres (receptor within 5 of ligand)")
cmd.show("sticks", "pocket and not (name C+N+O)")
cmd.color("gray70", "pocket and elem C")
cmd.util.cnc("pocket")
cmd.set("stick_radius", 0.13, "pocket")

# Ligand 1712 -- stable-donor BLUE carbons (#0072B2), heteroatoms by element.
cmd.show("sticks", "ligand")
cmd.set_color("stableblue", [0.0, 0.447, 0.698])
cmd.color("stableblue", "ligand and elem C")
cmd.util.cnc("ligand")
cmd.set("stick_radius", 0.22, "ligand")

# Framing: orient on the ligand, then pull back so nothing clips.
cmd.orient("ligand")
cmd.center("ligand or pocket")
cmd.zoom("ligand or pocket", 11)

# Quality
cmd.set("ray_opaque_background", 1)
cmd.set("antialias", 2)
cmd.set("ray_shadows", 0)
cmd.set("specular", 0.25)
cmd.ray(2400, 2000)
cmd.png(OUT_PNG, dpi=600)
print("wrote " + OUT_PNG)
