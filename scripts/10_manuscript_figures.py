#!/usr/bin/env python3
"""
Generate additional publication-quality figures for manuscript
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 11
plt.rcParams['font.family'] = 'serif'

ANALYSIS_DIR = Path('../analysis')
FIGURES_DIR = Path('../figures')
FIGURES_DIR.mkdir(exist_ok=True)

# Load data
reorg = pd.read_csv(ANALYSIS_DIR / 'reorganization_energies.csv')
ibo = pd.read_csv(ANALYSIS_DIR / 'ibo_analysis.csv')
stacking = pd.read_csv(ANALYSIS_DIR / 'stacking_analysis.csv')

reorg_ok = reorg[reorg['status'] == 'success'].copy()
ibo_ok = ibo[ibo['status'] == 'success'].copy()

print("Generating manuscript figures...")

# ============================================================================
# Figure: Reorganization Energies with Error Bars
# ============================================================================
fig, ax = plt.subplots(1, 1, figsize=(8, 6))

# Calculate statistics
lambda_hole_mean = reorg_ok['lambda_hole_eV'].mean()
lambda_hole_std = reorg_ok['lambda_hole_eV'].std()
lambda_e_mean = reorg_ok['lambda_electron_eV'].mean()
lambda_e_std = reorg_ok['lambda_electron_eV'].std()

x = np.arange(2)
means = [lambda_hole_mean, abs(lambda_e_mean)]
stds = [lambda_hole_std, lambda_e_std]
labels = ['λ$_{hole}$', 'λ$_{electron}$']
colors = ['#d62728', '#1f77b4']

bars = ax.bar(x, means, yerr=stds, capsize=10, color=colors, alpha=0.7, 
              edgecolor='black', linewidth=1.5)

ax.set_ylabel('Reorganization Energy (eV)', fontsize=12, fontweight='bold')
ax.set_title('Charge Transport Reorganization Energies', fontsize=13, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=12)
ax.grid(axis='y', alpha=0.3)

# Add value labels
for i, (m, s) in enumerate(zip(means, stds)):
    ax.text(i, m + s + 0.3, f'{m:.2f}±{s:.2f}', 
            ha='center', va='bottom', fontsize=10, fontweight='bold')

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'reorganization_energies_manuscript.png', dpi=300, bbox_inches='tight')
plt.savefig(FIGURES_DIR / 'reorganization_energies_manuscript.pdf', bbox_inches='tight')
print("  ✓ Saved: reorganization_energies_manuscript.png/pdf")
plt.close()

# ============================================================================
# Figure: IBO Localization Analysis
# ============================================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: Number of IBOs
axes[0].bar(range(len(ibo_ok)), ibo_ok['n_ibos'], color='#2ca02c', 
            alpha=0.7, edgecolor='black', linewidth=1.5)
axes[0].set_xlabel('Molecule ID', fontsize=12, fontweight='bold')
axes[0].set_ylabel('Number of IBOs', fontsize=12, fontweight='bold')
axes[0].set_title('(a) Intrinsic Bond Orbital Count', fontsize=12, fontweight='bold')
axes[0].set_xticks(range(len(ibo_ok)))
axes[0].set_xticklabels(ibo_ok['mol_id'], rotation=45)
axes[0].grid(axis='y', alpha=0.3)

# Add mean line
mean_ibos = ibo_ok['n_ibos'].mean()
axes[0].axhline(mean_ibos, color='red', linestyle='--', linewidth=2, 
                label=f'Mean: {mean_ibos:.1f}')
axes[0].legend()

# Panel B: Localization distance
axes[1].bar(range(len(ibo_ok)), ibo_ok['avg_ibo_atom_distance_angstrom'], 
            color='#ff7f0e', alpha=0.7, edgecolor='black', linewidth=1.5)
axes[1].set_xlabel('Molecule ID', fontsize=12, fontweight='bold')
axes[1].set_ylabel('IBO Localization Distance (Å)', fontsize=12, fontweight='bold')
axes[1].set_title('(b) IBO Spatial Localization', fontsize=12, fontweight='bold')
axes[1].set_xticks(range(len(ibo_ok)))
axes[1].set_xticklabels(ibo_ok['mol_id'], rotation=45)
axes[1].grid(axis='y', alpha=0.3)

# Add mean line
mean_loc = ibo_ok['avg_ibo_atom_distance_angstrom'].mean()
axes[1].axhline(mean_loc, color='red', linestyle='--', linewidth=2,
                label=f'Mean: {mean_loc:.2f} Å')
axes[1].legend()

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'ibo_analysis_manuscript.png', dpi=300, bbox_inches='tight')
plt.savefig(FIGURES_DIR / 'ibo_analysis_manuscript.pdf', bbox_inches='tight')
print("  ✓ Saved: ibo_analysis_manuscript.png/pdf")
plt.close()

# ============================================================================
# Figure: Combined Transport Properties
# ============================================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Panel A: Individual reorganization energies
mol_ids = reorg_ok['mol_id'].values
x_pos = np.arange(len(mol_ids))
width = 0.35

axes[0, 0].bar(x_pos - width/2, reorg_ok['lambda_hole_eV'], width, 
               label='λ$_{hole}$', color='#d62728', alpha=0.7, edgecolor='black')
axes[0, 0].bar(x_pos + width/2, abs(reorg_ok['lambda_electron_eV']), width,
               label='|λ$_{electron}$|', color='#1f77b4', alpha=0.7, edgecolor='black')
axes[0, 0].set_ylabel('Reorganization Energy (eV)', fontsize=11, fontweight='bold')
axes[0, 0].set_title('(a) Reorganization Energies by Molecule', fontsize=11, fontweight='bold')
axes[0, 0].set_xticks(x_pos)
axes[0, 0].set_xticklabels(mol_ids, rotation=45)
axes[0, 0].legend()
axes[0, 0].grid(axis='y', alpha=0.3)

# Panel B: Reorganization ratio
axes[0, 1].bar(range(len(reorg_ok)), abs(reorg_ok['lambda_ratio']), 
               color='#9467bd', alpha=0.7, edgecolor='black')
axes[0, 1].set_ylabel('|λ$_{hole}$/λ$_{electron}$|', fontsize=11, fontweight='bold')
axes[0, 1].set_title('(b) Hole/Electron Reorganization Ratio', fontsize=11, fontweight='bold')
axes[0, 1].set_xticks(range(len(reorg_ok)))
axes[0, 1].set_xticklabels(mol_ids, rotation=45)
axes[0, 1].set_yscale('log')
axes[0, 1].grid(axis='y', alpha=0.3)

# Panel C: IBO count vs localization
axes[1, 0].scatter(ibo_ok['n_ibos'], ibo_ok['avg_ibo_atom_distance_angstrom'],
                   s=150, c=range(len(ibo_ok)), cmap='viridis', 
                   alpha=0.7, edgecolor='black', linewidth=1.5)
axes[1, 0].set_xlabel('Number of IBOs', fontsize=11, fontweight='bold')
axes[1, 0].set_ylabel('Localization Distance (Å)', fontsize=11, fontweight='bold')
axes[1, 0].set_title('(c) IBO Count vs Localization', fontsize=11, fontweight='bold')
axes[1, 0].grid(alpha=0.3)

# Panel D: Summary statistics table
axes[1, 1].axis('off')
summary_text = f"""
Transport Properties Summary (n={len(reorg_ok)} molecules)

Reorganization Energies:
  λ_hole:     {lambda_hole_mean:.2f} ± {lambda_hole_std:.2f} eV
  λ_electron: {lambda_e_mean:.2f} ± {lambda_e_std:.2f} eV
  Ratio:      {abs(reorg_ok['lambda_ratio'].mean()):.1f} ± {reorg_ok['lambda_ratio'].std():.1f}

IBO Analysis (n={len(ibo_ok)} molecules):
  Average IBOs:    {ibo_ok['n_ibos'].mean():.1f} ± {ibo_ok['n_ibos'].std():.1f}
  Localization:    {mean_loc:.2f} ± {ibo_ok['avg_ibo_atom_distance_angstrom'].std():.2f} Å

Stacking Analysis (Mol 17851):
  Optimal distance: {stacking['optimal_separation_angstrom'].values[0]:.2f} Å
  Quality:          {stacking['stacking_quality'].values[0]}
"""
axes[1, 1].text(0.1, 0.5, summary_text, fontsize=10, family='monospace',
                verticalalignment='center', bbox=dict(boxstyle='round', 
                facecolor='wheat', alpha=0.3))
axes[1, 1].set_title('(d) Statistical Summary', fontsize=11, fontweight='bold')

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'transport_properties_combined.png', dpi=300, bbox_inches='tight')
plt.savefig(FIGURES_DIR / 'transport_properties_combined.pdf', bbox_inches='tight')
print("  ✓ Saved: transport_properties_combined.png/pdf")
plt.close()

print("\n" + "="*60)
print("All manuscript figures generated successfully!")
print("="*60)
print(f"\nFigures saved to: {FIGURES_DIR}")
print("\nGenerated files:")
print("  1. reorganization_energies_manuscript.png/pdf")
print("  2. ibo_analysis_manuscript.png/pdf")
print("  3. transport_properties_combined.png/pdf")
