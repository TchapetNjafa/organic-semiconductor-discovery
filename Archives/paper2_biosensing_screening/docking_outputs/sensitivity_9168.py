"""
Sensitivity analysis for Phase-5 review item C2.

9168 (C26H16) is a bare all-carbon polycyclic aromatic hydrocarbon: it has no
heteroatoms, unlike the other three "stable donors" (1712, 17574, 18506), which
are heteroatom-rich donor-acceptor structures. Yet 9168 is the single strongest
docking binder in the set (-8.0 kcal/mol) and thus anchors the "stable donors
bind strongly" conclusion. This script re-derives the three headline findings
WITH and WITHOUT 9168 to test whether the two-axis separation is an artifact of
including one structurally atypical molecule.

Outputs a compact markdown table to stdout and writes sensitivity_9168.csv.
"""
import pandas as pd

BASE = "/home/tchapet/Post-Doc/ARTICLES-DOSSIERS-DES-ARTICLES-EN-REDACTION/NOUVELS-AXES-DE-RECHERCHE-A-REGARDER-URGEMMENT/ARTICLES-EN-REDACTIONS/ARTICLE_MVOTO"
DIPOLE_CSV = f"{BASE}/MVOTOOPV/PUBCHEMQC_MULTIFONCTIONAL_2025/DATASET/PCE_paper_GDB9.csv"
DOCK_CSV = f"{BASE}/Paper2_JCIM/docking_rerun/docking_results_clean.csv"
SOLV_CAM = f"{BASE}/Paper2_JCIM/solvatochromic_tddft/solvatochromic_results_CAM-B3LYP.csv"
SOLV_B3 = f"{BASE}/Paper2_JCIM/solvatochromic_tddft/solvatochromic_results_B3LYP.csv"

STABLE = [1712, 17574, 18506, 9168]      # as published
REACTIVE = [20778, 17851, 4550]
OUTLIER = 9168


def s1_shift(csv):
    """Water-minus-toluene S1 (state 1) energy shift, eV, per molecule."""
    df = pd.read_csv(csv)
    s1 = df[df["state"] == 1]
    out = {}
    for m, g in s1.groupby("mol_id"):
        tol = g[g["solvent"] == "Toluene"]["energy_eV"].iloc[0]
        wat = g[g["solvent"] == "Water"]["energy_eV"].iloc[0]
        out[int(m)] = wat - tol
    return out


dip = pd.read_csv(DIPOLE_CSV).set_index("mol_id")["dipole_moment"].to_dict()
dock = pd.read_csv(DOCK_CSV)
best = dock.groupby("mol_id")["affinity_kcal_mol"].min().to_dict()   # strongest (most negative)
worst = dock.groupby("mol_id")["affinity_kcal_mol"].max().to_dict()  # weakest across targets
shift_cam = s1_shift(SOLV_CAM)
shift_b3 = s1_shift(SOLV_B3)

print("## Per-molecule summary (verify vs manuscript)\n")
print("| mol | class | dipole D | dock best | dock worst | |shift| CAM eV |")
print("|-----|-------|----------|-----------|------------|---------------|")
for m in STABLE + REACTIVE:
    cls = "stable" if m in STABLE else "reactive"
    tag = " (OUTLIER)" if m == OUTLIER else ""
    print(f"| {m}{tag} | {cls} | {dip[m]:.2f} | {best[m]:.1f} | {worst[m]:.1f} | {abs(shift_cam[m]):.3f} |")

def band(ids):
    b = [best[m] for m in ids]
    return min(b), max(b)  # strongest, weakest-of-strongest

print("\n## Finding 1 — docking separation by stability class\n")
for label, stable in [("WITH 9168", STABLE), ("WITHOUT 9168", [m for m in STABLE if m != OUTLIER])]:
    st_strong = min(best[m] for m in stable)   # best stable
    st_weak = max(best[m] for m in stable)     # weakest stable (boundary)
    rx_strong = min(best[m] for m in REACTIVE) # best reactive (boundary)
    rx_weak = max(best[m] for m in REACTIVE)
    gap = rx_strong - st_weak                  # >0 means clean separation
    print(f"**{label}**: stable best-affinity band [{st_strong:.1f}, {st_weak:.1f}] | "
          f"reactive [{rx_strong:.1f}, {rx_weak:.1f}] | "
          f"boundary gap (reactive-best minus stable-weakest) = {gap:+.1f} kcal/mol "
          f"-> {'clean separation' if gap > 0 else 'OVERLAP of ' + format(-gap,'.1f')}")

# per-target check without 9168
print("\n### Per-target boundary (WITHOUT 9168): weakest stable vs strongest reactive\n")
print("| target | weakest stable | strongest reactive | separated? |")
print("|--------|----------------|--------------------|-----------|")
for tgt, g in dock.groupby("target"):
    gg = g.set_index("mol_id")["affinity_kcal_mol"].to_dict()
    st = [gg[m] for m in STABLE if m != OUTLIER and m in gg]
    rx = [gg[m] for m in REACTIVE if m in gg]
    ws, sr = max(st), min(rx)
    print(f"| {tgt} | {ws:.1f} | {sr:.1f} | {'yes' if sr > ws else 'NO ('+format(ws-sr,'.1f')+' overlap)'} |")

print("\n## Finding 2 — dipole partition\n")
for label, stable in [("WITH 9168", STABLE), ("WITHOUT 9168", [m for m in STABLE if m != OUTLIER])]:
    st_d = [dip[m] for m in stable]
    rx_d = [dip[m] for m in REACTIVE]
    print(f"**{label}**: stable dipole [{min(st_d):.2f}, {max(st_d):.2f}] D | "
          f"reactive [{min(rx_d):.2f}, {max(rx_d):.2f}] D | "
          f"gap = {min(rx_d) - max(st_d):+.2f} D")

print("\n## Finding 3 — solvatochromic shift (CAM-B3LYP |shift|)\n")
for label, stable in [("WITH 9168", STABLE), ("WITHOUT 9168", [m for m in STABLE if m != OUTLIER])]:
    st_s = [abs(shift_cam[m]) for m in stable]
    rx_s = [abs(shift_cam[m]) for m in REACTIVE]
    print(f"**{label}**: stable |shift| max = {max(st_s):.3f} eV | "
          f"reactive |shift| min = {min(rx_s):.3f} eV | "
          f"gap = {min(rx_s) - max(st_s):+.3f} eV")

# write csv
rows = []
for m in STABLE + REACTIVE:
    rows.append({"mol_id": m, "class": "stable" if m in STABLE else "reactive",
                 "is_9168_outlier": m == OUTLIER, "dipole_D": round(dip[m], 3),
                 "dock_best_kcal": best[m], "dock_worst_kcal": worst[m],
                 "shift_CAM_eV": round(shift_cam[m], 4), "shift_B3_eV": round(shift_b3[m], 4)})
pd.DataFrame(rows).to_csv(f"{BASE}/Paper2_JCIM/docking_rerun/sensitivity_9168.csv", index=False)
print("\nWrote sensitivity_9168.csv")
