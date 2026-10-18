#!/usr/bin/env python3
"""
Comprehensive Analysis and Visualization of NTO Study Results
Generates publication-quality figures and statistical summaries
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 10

# Directories
ANALYSIS_DIR = Path('../analysis')
FIGURES_DIR = Path('../figures')
FIGURES_DIR.mkdir(exist_ok=True)

# Load all datasets
print("Loading datasets...")
xtb = pd.read_csv(ANALYSIS_DIR / 'xtb_excited_state_screening.csv')
multiwfn = pd.read_csv(ANALYSIS_DIR / 'multiwfn_nto_analysis.csv')
pyscf = pd.read_csv(ANALYSIS_DIR / 'pyscf_detailed_analysis.csv')
reorg = pd.read_csv(ANALYSIS_DIR / 'reorganization_energies.csv')
stacking = pd.read_csv(ANALYSIS_DIR / 'stacking_analysis.csv')
nbo = pd.read_csv(ANALYSIS_DIR / 'nbo_analysis.csv')
ibo = pd.read_csv(ANALYSIS_DIR / 'ibo_analysis.csv')

print(f"Loaded: xTB ({len(xtb)}), Multiwfn ({len(multiwfn)}), PySCF ({len(pyscf)})")
print(f"        Reorg ({len(reorg)}), Stacking ({len(stacking)}), NBO ({len(nbo)}), IBO ({len(ibo)})")

# Filter successful calculations
xtb_ok = xtb[xtb['status'] == 'success'].copy()
multiwfn_ok = multiwfn[multiwfn['status'] == 'success'].copy()
pyscf_ok = pyscf[pyscf['status'] == 'success'].copy()
reorg_ok = reorg[reorg['status'] == 'success'].copy()

print(f"\nSuccessful: xTB ({len(xtb_ok)}), Multiwfn ({len(multiwfn_ok)}), PySCF ({len(pyscf_ok)}), Reorg ({len(reorg_ok)})")

# ============================================================================
# Figure 1: xTB Excited State Screening Overview
# ============================================================================
print("\nGenerating Figure 1: xTB Screening Overview...")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# S1 energy distribution
axes[0, 0].hist(xtb_ok['S1_eV'], bins=50, color='steelblue', alpha=0.7, edgecolor='black')
axes[0, 0].axvline(1.5, color='red', linestyle='--', label='OPV range')
axes[0, 0].axvline(3.0, color='red', linestyle='--')
axes[0, 0].set_xlabel('S1 Energy (eV)')
axes[0, 0].set_ylabel('Count')
axes[0, 0].set_title(f'S1 Excitation Energy Distribution (n={len(xtb_ok)})')
axes[0, 0].legend()

# Oscillator strength
axes[0, 1].hist(xtb_ok['S1_f'], bins=50, color='coral', alpha=0.7, edgecolor='black')
axes[0, 1].axvline(0.1, color='red', linestyle='--', label='Strong absorption')
axes[0, 1].set_xlabel('S1 Oscillator Strength')
axes[0, 1].set_ylabel('Count')
axes[0, 1].set_title('S1 Oscillator Strength Distribution')
axes[0, 1].legend()

# S1 vs S2 energy
axes[1, 0].scatter(xtb_ok['S1_eV'], xtb_ok['S2_eV'], alpha=0.3, s=10, c='purple')
axes[1, 0].plot([0, 10], [0, 10], 'k--', alpha=0.3)
axes[1, 0].set_xlabel('S1 Energy (eV)')
axes[1, 0].set_ylabel('S2 Energy (eV)')
axes[1, 0].set_title('S1 vs S2 Excitation Energies')

# S1 energy vs oscillator strength
scatter = axes[1, 1].scatter(xtb_ok['S1_eV'], xtb_ok['S1_f'], 
                             c=xtb_ok['S1_nm'], cmap='viridis', alpha=0.5, s=10)
axes[1, 1].set_xlabel('S1 Energy (eV)')
axes[1, 1].set_ylabel('S1 Oscillator Strength')
axes[1, 1].set_title('S1 Energy vs Oscillator Strength')
cbar = plt.colorbar(scatter, ax=axes[1, 1])
cbar.set_label('Wavelength (nm)')

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'fig1_xtb_screening_overview.png', dpi=300, bbox_inches='tight')
print(f"  Saved: fig1_xtb_screening_overview.png")
plt.close()

# ============================================================================
# Figure 2: Multiwfn NTO Analysis
# ============================================================================
print("\nGenerating Figure 2: NTO Analysis...")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Lambda_1 distribution
if 'Lambda_1' in multiwfn_ok.columns:
    axes[0, 0].hist(multiwfn_ok['Lambda_1'].dropna(), bins=50, color='green', alpha=0.7, edgecolor='black')
    axes[0, 0].axvline(0.8, color='red', linestyle='--', label='Single excitation')
    axes[0, 0].set_xlabel('Lambda_1 (NTO weight)')
    axes[0, 0].set_ylabel('Count')
    axes[0, 0].set_title('NTO Lambda_1 Distribution')
    axes[0, 0].legend()

# Sr index (electron-hole overlap)
if 'Sr_index' in multiwfn_ok.columns:
    axes[0, 1].hist(multiwfn_ok['Sr_index'].dropna(), bins=50, color='orange', alpha=0.7, edgecolor='black')
    axes[0, 1].set_xlabel('Sr Index (e-h overlap)')
    axes[0, 1].set_ylabel('Count')
    axes[0, 1].set_title('Electron-Hole Overlap Distribution')

# D index (CT distance)
if 'D_index_Angstrom' in multiwfn_ok.columns:
    axes[1, 0].hist(multiwfn_ok['D_index_Angstrom'].dropna(), bins=50, color='teal', alpha=0.7, edgecolor='black')
    axes[1, 0].axvline(2.0, color='red', linestyle='--', label='CT character')
    axes[1, 0].set_xlabel('D Index (Å)')
    axes[1, 0].set_ylabel('Count')
    axes[1, 0].set_title('Charge Transfer Distance Distribution')
    axes[1, 0].legend()

# Sr vs D index
if 'Sr_index' in multiwfn_ok.columns and 'D_index_Angstrom' in multiwfn_ok.columns:
    axes[1, 1].scatter(multiwfn_ok['D_index_Angstrom'], multiwfn_ok['Sr_index'], 
                       alpha=0.3, s=10, c='purple')
    axes[1, 1].set_xlabel('D Index (Å)')
    axes[1, 1].set_ylabel('Sr Index')
    axes[1, 1].set_title('CT Distance vs e-h Overlap')

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'fig2_nto_analysis.png', dpi=300, bbox_inches='tight')
print(f"  Saved: fig2_nto_analysis.png")
plt.close()

# ============================================================================
# Figure 3: Top Molecules Comparison
# ============================================================================
print("\nGenerating Figure 3: Top Molecules Comparison...")

# Merge datasets for top molecules
top_mols = [17851, 19531, 14380, 20778, 977, 7801, 18985]
top_data = []

for mol_id in top_mols:
    row = {'mol_id': mol_id}
    
    # xTB data
    xtb_row = xtb_ok[xtb_ok['mol_id'] == mol_id]
    if len(xtb_row) > 0:
        row['S1_eV'] = xtb_row['S1_eV'].values[0]
        row['S1_f'] = xtb_row['S1_f'].values[0]
    
    # Reorganization
    reorg_row = reorg_ok[reorg_ok['mol_id'] == mol_id]
    if len(reorg_row) > 0:
        row['lambda_hole'] = reorg_row['lambda_hole_eV'].values[0]
        row['lambda_electron'] = reorg_row['lambda_electron_eV'].values[0]
    
    # IBO
    ibo_row = ibo[ibo['mol_id'] == mol_id]
    if len(ibo_row) > 0:
        row['n_ibos'] = ibo_row['n_ibos'].values[0]
        row['ibo_localization'] = ibo_row['avg_ibo_atom_distance_angstrom'].values[0]
    
    top_data.append(row)

top_df = pd.DataFrame(top_data)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# S1 energies
if 'S1_eV' in top_df.columns:
    axes[0, 0].bar(range(len(top_df)), top_df['S1_eV'], color='steelblue', edgecolor='black')
    axes[0, 0].set_xticks(range(len(top_df)))
    axes[0, 0].set_xticklabels(top_df['mol_id'], rotation=45)
    axes[0, 0].set_ylabel('S1 Energy (eV)')
    axes[0, 0].set_title('S1 Excitation Energies')
    axes[0, 0].grid(axis='y', alpha=0.3)

# Oscillator strengths
if 'S1_f' in top_df.columns:
    axes[0, 1].bar(range(len(top_df)), top_df['S1_f'], color='coral', edgecolor='black')
    axes[0, 1].set_xticks(range(len(top_df)))
    axes[0, 1].set_xticklabels(top_df['mol_id'], rotation=45)
    axes[0, 1].set_ylabel('Oscillator Strength')
    axes[0, 1].set_title('S1 Oscillator Strengths')
    axes[0, 1].grid(axis='y', alpha=0.3)

# Reorganization energies
if 'lambda_hole' in top_df.columns and 'lambda_electron' in top_df.columns:
    x = np.arange(len(top_df))
    width = 0.35
    axes[1, 0].bar(x - width/2, top_df['lambda_hole'], width, label='λ_hole', color='red', alpha=0.7)
    axes[1, 0].bar(x + width/2, top_df['lambda_electron'], width, label='λ_electron', color='blue', alpha=0.7)
    axes[1, 0].set_xticks(x)
    axes[1, 0].set_xticklabels(top_df['mol_id'], rotation=45)
    axes[1, 0].set_ylabel('Reorganization Energy (eV)')
    axes[1, 0].set_title('Reorganization Energies')
    axes[1, 0].legend()
    axes[1, 0].set_yscale('log')
    axes[1, 0].grid(axis='y', alpha=0.3)

# IBO analysis
if 'n_ibos' in top_df.columns:
    axes[1, 1].bar(range(len(top_df)), top_df['n_ibos'], color='green', edgecolor='black')
    axes[1, 1].set_xticks(range(len(top_df)))
    axes[1, 1].set_xticklabels(top_df['mol_id'], rotation=45)
    axes[1, 1].set_ylabel('Number of IBOs')
    axes[1, 1].set_title('Intrinsic Bond Orbitals')
    axes[1, 1].grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'fig3_top_molecules_comparison.png', dpi=300, bbox_inches='tight')
print(f"  Saved: fig3_top_molecules_comparison.png")
plt.close()

# ============================================================================
# Summary Statistics
# ============================================================================
print("\n" + "="*60)
print("SUMMARY STATISTICS")
print("="*60)

print("\n1. xTB Screening (17,458 molecules):")
print(f"   Success rate: {len(xtb_ok)/len(xtb)*100:.1f}%")
print(f"   S1 energy: {xtb_ok['S1_eV'].mean():.2f} ± {xtb_ok['S1_eV'].std():.2f} eV")
print(f"   S1 f: {xtb_ok['S1_f'].mean():.3f} ± {xtb_ok['S1_f'].std():.3f}")
print(f"   OPV range (1.5-3.0 eV): {((xtb_ok['S1_eV'] >= 1.5) & (xtb_ok['S1_eV'] <= 3.0)).sum()} molecules")

print("\n2. Multiwfn NTO Analysis:")
print(f"   Success rate: {len(multiwfn_ok)/len(multiwfn)*100:.1f}%")
if 'Lambda_1' in multiwfn_ok.columns:
    print(f"   Lambda_1: {multiwfn_ok['Lambda_1'].mean():.3f} ± {multiwfn_ok['Lambda_1'].std():.3f}")
if 'Sr_index' in multiwfn_ok.columns:
    print(f"   Sr index: {multiwfn_ok['Sr_index'].mean():.3f} ± {multiwfn_ok['Sr_index'].std():.3f}")
if 'D_index_Angstrom' in multiwfn_ok.columns:
    print(f"   D index: {multiwfn_ok['D_index_Angstrom'].mean():.2f} ± {multiwfn_ok['D_index_Angstrom'].std():.2f} Å")

print("\n3. Reorganization Energies (7 molecules):")
print(f"   λ_hole: {reorg_ok['lambda_hole_eV'].mean():.2f} ± {reorg_ok['lambda_hole_eV'].std():.2f} eV")
print(f"   λ_electron: {reorg_ok['lambda_electron_eV'].mean():.2f} ± {reorg_ok['lambda_electron_eV'].std():.2f} eV")
print(f"   λ_ratio: {reorg_ok['lambda_ratio'].mean():.1f} ± {reorg_ok['lambda_ratio'].std():.1f}")

print("\n4. Stacking Analysis:")
print(f"   Molecule 17851: {stacking['optimal_separation_angstrom'].values[0]:.2f} Å, {stacking['stacking_quality'].values[0]}")

print("\n5. IBO Analysis:")
print(f"   Molecules analyzed: {len(ibo)}")
print(f"   Average IBOs: {ibo['n_ibos'].mean():.1f}")
print(f"   Localization: {ibo['avg_ibo_atom_distance_angstrom'].mean():.3f} Å")

print("\n" + "="*60)
print("Analysis complete! Figures saved to:", FIGURES_DIR)
print("="*60)
