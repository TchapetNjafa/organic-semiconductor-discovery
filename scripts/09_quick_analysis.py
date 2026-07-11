#!/usr/bin/env python3
"""
Quick Analysis of Available NTO Study Results
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300

ANALYSIS_DIR = Path('../analysis')
FIGURES_DIR = Path('../figures')
FIGURES_DIR.mkdir(exist_ok=True)

# Load datasets
print("Loading datasets...")
pyscf = pd.read_csv(ANALYSIS_DIR / 'pyscf_detailed_analysis.csv')
reorg = pd.read_csv(ANALYSIS_DIR / 'reorganization_energies.csv')
stacking = pd.read_csv(ANALYSIS_DIR / 'stacking_analysis.csv')
nbo = pd.read_csv(ANALYSIS_DIR / 'nbo_analysis.csv')
ibo = pd.read_csv(ANALYSIS_DIR / 'ibo_analysis.csv')

pyscf_ok = pyscf[pyscf['status'] == 'success'].copy()
reorg_ok = reorg[reorg['status'] == 'success'].copy()

print(f"PySCF: {len(pyscf_ok)}/{len(pyscf)} successful")
print(f"Reorg: {len(reorg_ok)}/{len(reorg)} successful")
print(f"Stacking: {len(stacking)} molecules")
print(f"NBO: {len(nbo)} molecules")
print(f"IBO: {len(ibo)} molecules")

# Create summary table
top_mols = [17851, 19531, 14380, 20778, 977, 7801, 18985]
summary = []

for mol_id in top_mols:
    row = {'Molecule': mol_id}
    
    # PySCF
    p = pyscf_ok[pyscf_ok['mol_id'] == mol_id]
    if len(p) > 0:
        row['S1 (eV)'] = f"{p['S1_eV'].values[0]:.2f}"
        row['S1 f'] = f"{p['S1_f'].values[0]:.3f}"
        row['HOMO (eV)'] = f"{p['E_HOMO'].values[0]:.2f}"
        row['LUMO (eV)'] = f"{p['E_LUMO'].values[0]:.2f}"
    
    # Reorganization
    r = reorg_ok[reorg_ok['mol_id'] == mol_id]
    if len(r) > 0:
        row['λ_hole (eV)'] = f"{r['lambda_hole_eV'].values[0]:.2f}"
        row['λ_e (eV)'] = f"{r['lambda_electron_eV'].values[0]:.2f}"
    
    # IBO
    i = ibo[ibo['mol_id'] == mol_id]
    if len(i) > 0 and not pd.isna(i['n_ibos'].values[0]):
        row['IBOs'] = int(i['n_ibos'].values[0])
        row['IBO_loc (Å)'] = f"{i['avg_ibo_atom_distance_angstrom'].values[0]:.2f}"
    
    summary.append(row)

summary_df = pd.DataFrame(summary)
print("\n" + "="*80)
print("TOP MOLECULES SUMMARY")
print("="*80)
print(summary_df.to_string(index=False))

# Save summary
summary_df.to_csv(ANALYSIS_DIR / 'top_molecules_summary.csv', index=False)
print(f"\nSaved: {ANALYSIS_DIR / 'top_molecules_summary.csv'}")

# Create figure
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Reorganization energies
if len(reorg_ok) > 0:
    x = np.arange(len(reorg_ok))
    axes[0, 0].bar(x - 0.2, reorg_ok['lambda_hole_eV'], 0.4, label='λ_hole', color='red', alpha=0.7)
    axes[0, 0].bar(x + 0.2, reorg_ok['lambda_electron_eV'], 0.4, label='λ_electron', color='blue', alpha=0.7)
    axes[0, 0].set_xticks(x)
    axes[0, 0].set_xticklabels(reorg_ok['mol_id'], rotation=45)
    axes[0, 0].set_ylabel('Reorganization Energy (eV)')
    axes[0, 0].set_title('Reorganization Energies')
    axes[0, 0].legend()
    axes[0, 0].set_yscale('log')

# IBO counts
if len(ibo) > 0:
    axes[0, 1].bar(range(len(ibo)), ibo['n_ibos'], color='green', edgecolor='black')
    axes[0, 1].set_xticks(range(len(ibo)))
    axes[0, 1].set_xticklabels(ibo['mol_id'], rotation=45)
    axes[0, 1].set_ylabel('Number of IBOs')
    axes[0, 1].set_title('Intrinsic Bond Orbitals')

# PySCF S1 energies
if len(pyscf_ok) > 0:
    axes[1, 0].bar(range(len(pyscf_ok)), pyscf_ok['S1_eV'], color='steelblue', edgecolor='black')
    axes[1, 0].set_xticks(range(len(pyscf_ok)))
    axes[1, 0].set_xticklabels(pyscf_ok['mol_id'], rotation=45)
    axes[1, 0].set_ylabel('S1 Energy (eV)')
    axes[1, 0].set_title('S1 Excitation Energies (CAM-B3LYP)')

# PySCF oscillator strengths
if len(pyscf_ok) > 0:
    axes[1, 1].bar(range(len(pyscf_ok)), pyscf_ok['S1_f'], color='coral', edgecolor='black')
    axes[1, 1].set_xticks(range(len(pyscf_ok)))
    axes[1, 1].set_xticklabels(pyscf_ok['mol_id'], rotation=45)
    axes[1, 1].set_ylabel('Oscillator Strength')
    axes[1, 1].set_title('S1 Oscillator Strengths')

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'top_molecules_analysis.png', dpi=300, bbox_inches='tight')
print(f"Saved: {FIGURES_DIR / 'top_molecules_analysis.png'}")

print("\n" + "="*80)
print("STATISTICS")
print("="*80)
print(f"\nReorganization Energies (n={len(reorg_ok)}):")
print(f"  λ_hole: {reorg_ok['lambda_hole_eV'].mean():.2f} ± {reorg_ok['lambda_hole_eV'].std():.2f} eV")
print(f"  λ_electron: {reorg_ok['lambda_electron_eV'].mean():.2f} ± {reorg_ok['lambda_electron_eV'].std():.2f} eV")
print(f"  Ratio: {reorg_ok['lambda_ratio'].mean():.1f} ± {reorg_ok['lambda_ratio'].std():.1f}")

print(f"\nIBO Analysis (n={len(ibo)}):")
print(f"  Average IBOs: {ibo['n_ibos'].mean():.1f} ± {ibo['n_ibos'].std():.1f}")
print(f"  Localization: {ibo['avg_ibo_atom_distance_angstrom'].mean():.2f} ± {ibo['avg_ibo_atom_distance_angstrom'].std():.2f} Å")

print(f"\nStacking Analysis:")
print(f"  Mol 17851: {stacking['optimal_separation_angstrom'].values[0]:.2f} Å ({stacking['stacking_quality'].values[0]})")

print("\nAnalysis complete!")
