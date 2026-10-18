#!/usr/bin/env python3
"""
Regenerate all Phase 12-14 figures with expanded data (n=13-17)

MANDATORY REQUIREMENTS (from MVOTO_ARTICLE_ENHANCEMENT.md):
1. Use scienceplots package
2. Generate 3 formats: PNG (300 dpi), PDF, PGF
3. Label subplots: (a), (b), (c), etc.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend (REQUIRED for server)
import matplotlib.pyplot as plt
import scienceplots
from pathlib import Path

# Use scienceplots style (MANDATORY)
plt.style.use(['science', 'no-latex'])

BASE_DIR = Path(__file__).parent.parent.resolve()
DATA_DIR = BASE_DIR / 'analysis'
FIG_DIR = BASE_DIR / 'figures'
FIG_DIR.mkdir(exist_ok=True)

def save_figure(fig, name, base_dir=FIG_DIR):
    """Save figure in 3 required formats (MANDATORY)"""
    for fmt in ['png', 'pdf', 'pgf']:
        output_path = base_dir / f'{name}.{fmt}'
        if fmt == 'png':
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
        else:
            fig.savefig(output_path, bbox_inches='tight')
    print(f"✓ {name} saved (PNG, PDF, PGF)")

def add_subplot_label(ax, label, x=-0.15, y=1.05):
    """Add subplot label (a), (b), (c) - MANDATORY for multi-panel figures"""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=12, fontweight='bold')

# Load expanded data
print("Loading expanded data...")
orbital = pd.read_csv(DATA_DIR / 'orbital_composition_expanded.csv')
reorg = pd.read_csv(DATA_DIR / 'reorganization_energies_expanded.csv')
stack = pd.read_csv(DATA_DIR / 'stacking_analysis_expanded.csv')

print(f"  Orbital: n={len(orbital)}")
print(f"  Reorganization: n={len(reorg)}")
print(f"  Stacking: n={len(stack)}")

# ============================================================================
# FIGURE 1: Reorganization Energies (n=13)
# ============================================================================
print("\nGenerating Figure 1: Reorganization energies...")

fig, ax = plt.subplots(1, 1, figsize=(8, 5))

# Sort by lambda_hole for better visualization
reorg_sorted = reorg.sort_values('lambda_hole_eV')
mol_ids = reorg_sorted['mol_id'].astype(str)
lambdas = reorg_sorted['lambda_hole_eV']

# Color code: green if < 0.5 eV (excellent), orange otherwise
colors = ['#2ecc71' if l < 0.5 else '#e67e22' for l in lambdas]

bars = ax.bar(range(len(mol_ids)), lambdas, color=colors, alpha=0.8, edgecolor='black')

# Add threshold line
ax.axhline(y=0.5, color='red', linestyle='--', linewidth=1.5, 
           label='Excellent transport threshold')

ax.set_xlabel('Molecule ID', fontsize=11)
ax.set_ylabel('λ$_{hole}$ (eV)', fontsize=11)
ax.set_title('Reorganization Energy for Hole Transport (n=13)', fontsize=12, fontweight='bold')
ax.set_xticks(range(len(mol_ids)))
ax.set_xticklabels(mol_ids, rotation=45, ha='right')
ax.legend(loc='upper left', frameon=True)
ax.grid(axis='y', alpha=0.3)

# Add statistics text
stats_text = f"Mean: {lambdas.mean():.3f} ± {lambdas.std():.3f} eV\n"
stats_text += f"Range: {lambdas.min():.3f} - {lambdas.max():.3f} eV\n"
stats_text += f"λ < 0.5 eV: {(lambdas < 0.5).sum()}/{len(lambdas)}"
ax.text(0.98, 0.97, stats_text, transform=ax.transAxes,
        fontsize=9, verticalalignment='top', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

plt.tight_layout()
save_figure(fig, 'reorganization_energies_expanded')
plt.close()

# ============================================================================
# FIGURE 2: Orbital Composition (n=17)
# ============================================================================
print("Generating Figure 2: Orbital composition...")

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Sort by HOMO_N for better visualization
orbital_sorted = orbital.sort_values('HOMO_N_percent', ascending=False)
mol_ids = orbital_sorted['mol_id'].astype(str)

# Panel (a): HOMO composition
ax = axes[0]
homo_n = orbital_sorted['HOMO_N_percent']
homo_c = orbital_sorted['HOMO_C_percent']
homo_o = orbital_sorted['HOMO_O_percent']
homo_other = 100 - (homo_n + homo_c + homo_o)

bottom = np.zeros(len(mol_ids))
ax.bar(range(len(mol_ids)), homo_n, label='N', color='#3498db', alpha=0.8, bottom=bottom)
bottom += homo_n
ax.bar(range(len(mol_ids)), homo_c, label='C', color='#95a5a6', alpha=0.8, bottom=bottom)
bottom += homo_c
ax.bar(range(len(mol_ids)), homo_o, label='O', color='#e74c3c', alpha=0.8, bottom=bottom)
bottom += homo_o
ax.bar(range(len(mol_ids)), homo_other, label='Others', color='#f39c12', alpha=0.8, bottom=bottom)

ax.set_xlabel('Molecule ID', fontsize=11)
ax.set_ylabel('Composition (%)', fontsize=11)
ax.set_title('HOMO Composition (n=17)', fontsize=12, fontweight='bold')
ax.set_xticks(range(len(mol_ids)))
ax.set_xticklabels(mol_ids, rotation=45, ha='right', fontsize=8)
ax.legend(loc='upper right', frameon=True)
ax.grid(axis='y', alpha=0.3)
add_subplot_label(ax, '(a)')

# Panel (b): LUMO composition
ax = axes[1]
lumo_n = orbital_sorted['LUMO_N_percent']
lumo_c = orbital_sorted['LUMO_C_percent']
lumo_o = orbital_sorted['LUMO_O_percent']
lumo_other = 100 - (lumo_n + lumo_c + lumo_o)

bottom = np.zeros(len(mol_ids))
ax.bar(range(len(mol_ids)), lumo_n, label='N', color='#3498db', alpha=0.8, bottom=bottom)
bottom += lumo_n
ax.bar(range(len(mol_ids)), lumo_c, label='C', color='#95a5a6', alpha=0.8, bottom=bottom)
bottom += lumo_c
ax.bar(range(len(mol_ids)), lumo_o, label='O', color='#e74c3c', alpha=0.8, bottom=bottom)
bottom += lumo_o
ax.bar(range(len(mol_ids)), lumo_other, label='Others', color='#f39c12', alpha=0.8, bottom=bottom)

ax.set_xlabel('Molecule ID', fontsize=11)
ax.set_ylabel('Composition (%)', fontsize=11)
ax.set_title('LUMO Composition (n=17)', fontsize=12, fontweight='bold')
ax.set_xticks(range(len(mol_ids)))
ax.set_xticklabels(mol_ids, rotation=45, ha='right', fontsize=8)
ax.legend(loc='upper right', frameon=True)
ax.grid(axis='y', alpha=0.3)
add_subplot_label(ax, '(b)')

plt.tight_layout()
save_figure(fig, 'orbital_composition_expanded')
plt.close()

# ============================================================================
# FIGURE 3: Stacking Energy Curves (n=13)
# ============================================================================
print("Generating Figure 3: Stacking energy curves...")

fig, ax = plt.subplots(1, 1, figsize=(10, 6))

# Parse distance and energy lists
stack_sorted = stack.sort_values('interaction_energy_kcal_mol')

# Select representative molecules: top 3 attractive, 3 repulsive, 2 neutral
n_attractive = min(3, (stack['interaction_energy_kcal_mol'] < -0.5).sum())
n_repulsive = 3
attractive = stack_sorted.head(n_attractive)
repulsive = stack_sorted.tail(n_repulsive)

import ast
for idx, row in attractive.iterrows():
    distances = ast.literal_eval(row['distances'])
    energies = ast.literal_eval(row['energies'])
    ax.plot(distances, energies, 'o-', linewidth=2, markersize=6,
            label=f"{row['mol_id']} (attractive)")

for idx, row in repulsive.iterrows():
    distances = ast.literal_eval(row['distances'])
    energies = ast.literal_eval(row['energies'])
    ax.plot(distances, energies, 's--', linewidth=1.5, markersize=5, alpha=0.7,
            label=f"{row['mol_id']} (repulsive)")

ax.axhline(y=0, color='black', linestyle='-', linewidth=1, alpha=0.5)
ax.set_xlabel('Intermolecular Distance (Å)', fontsize=11)
ax.set_ylabel('Interaction Energy (kcal/mol)', fontsize=11)
ax.set_title('Stacking Interaction Energy Curves (n=13)', fontsize=12, fontweight='bold')
ax.legend(loc='best', frameon=True, fontsize=8, ncol=2)
ax.grid(alpha=0.3)

# Add statistics
stats_text = f"Mean: {stack['interaction_energy_kcal_mol'].mean():.2f} ± "
stats_text += f"{stack['interaction_energy_kcal_mol'].std():.2f} kcal/mol\n"
stats_text += f"Attractive (<0): {(stack['interaction_energy_kcal_mol'] < 0).sum()}/{len(stack)}"
ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
        fontsize=9, verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

plt.tight_layout()
save_figure(fig, 'stacking_energy_curves_expanded')
plt.close()

# ============================================================================
# FIGURE 4: Property Correlations (n=13 merged)
# ============================================================================
print("Generating Figure 4: Property correlations...")

# Merge datasets
merged = reorg.merge(orbital, on='mol_id')
merged = merged.merge(stack, on='mol_id')

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Panel (a): N in HOMO vs λ_hole
ax = axes[0]
ax.scatter(merged['HOMO_N_percent'], merged['lambda_hole_eV'], 
           s=100, alpha=0.7, edgecolor='black')
for idx, row in merged.iterrows():
    ax.annotate(str(row['mol_id']), 
                (row['HOMO_N_percent'], row['lambda_hole_eV']),
                fontsize=7, alpha=0.7)

from scipy.stats import pearsonr
r, p = pearsonr(merged['HOMO_N_percent'], merged['lambda_hole_eV'])
ax.set_xlabel('N in HOMO (%)', fontsize=11)
ax.set_ylabel('λ$_{hole}$ (eV)', fontsize=11)
ax.set_title(f'N in HOMO vs λ$_{{hole}}$ (n={len(merged)})', fontsize=11, fontweight='bold')
ax.text(0.05, 0.95, f'r = {r:.3f}\np = {p:.3f}', transform=ax.transAxes,
        fontsize=9, verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
ax.grid(alpha=0.3)
add_subplot_label(ax, '(a)')

# Panel (b): N in LUMO vs λ_hole
ax = axes[1]
ax.scatter(merged['LUMO_N_percent'], merged['lambda_hole_eV'],
           s=100, alpha=0.7, edgecolor='black', color='orange')
for idx, row in merged.iterrows():
    ax.annotate(str(row['mol_id']),
                (row['LUMO_N_percent'], row['lambda_hole_eV']),
                fontsize=7, alpha=0.7)

r, p = pearsonr(merged['LUMO_N_percent'], merged['lambda_hole_eV'])
ax.set_xlabel('N in LUMO (%)', fontsize=11)
ax.set_ylabel('λ$_{hole}$ (eV)', fontsize=11)
ax.set_title(f'N in LUMO vs λ$_{{hole}}$ (n={len(merged)})', fontsize=11, fontweight='bold')
ax.text(0.05, 0.95, f'r = {r:.3f}\np = {p:.3f}', transform=ax.transAxes,
        fontsize=9, verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
ax.grid(alpha=0.3)
add_subplot_label(ax, '(b)')

# Panel (c): λ_hole vs stacking energy
ax = axes[2]
ax.scatter(merged['lambda_hole_eV'], merged['interaction_energy_kcal_mol'],
           s=100, alpha=0.7, edgecolor='black', color='green')
for idx, row in merged.iterrows():
    ax.annotate(str(row['mol_id']),
                (row['lambda_hole_eV'], row['interaction_energy_kcal_mol']),
                fontsize=7, alpha=0.7)

r, p = pearsonr(merged['lambda_hole_eV'], merged['interaction_energy_kcal_mol'])
ax.set_xlabel('λ$_{hole}$ (eV)', fontsize=11)
ax.set_ylabel('Stacking Energy (kcal/mol)', fontsize=11)
ax.set_title(f'λ$_{{hole}}$ vs Stacking (n={len(merged)})', fontsize=11, fontweight='bold')
ax.text(0.05, 0.95, f'r = {r:.3f}\np = {p:.3f}', transform=ax.transAxes,
        fontsize=9, verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
ax.axhline(y=0, color='red', linestyle='--', linewidth=1, alpha=0.5)
ax.grid(alpha=0.3)
add_subplot_label(ax, '(c)')

plt.tight_layout()
save_figure(fig, 'property_correlations_expanded')
plt.close()

# ============================================================================
# FIGURE 5: Statistical Summary (n=13-17)
# ============================================================================
print("Generating Figure 5: Statistical summary...")

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Panel (a): Reorganization energy distribution
ax = axes[0]
ax.hist(reorg['lambda_hole_eV'], bins=8, color='#3498db', alpha=0.7, edgecolor='black')
ax.axvline(x=reorg['lambda_hole_eV'].mean(), color='red', linestyle='--', 
           linewidth=2, label='Mean')
ax.axvline(x=0.5, color='orange', linestyle='--', linewidth=2, 
           label='Threshold (0.5 eV)')
ax.set_xlabel('λ$_{hole}$ (eV)', fontsize=11)
ax.set_ylabel('Frequency', fontsize=11)
ax.set_title(f'Reorganization Energy Distribution (n={len(reorg)})', 
             fontsize=11, fontweight='bold')
ax.legend(frameon=True)
ax.grid(axis='y', alpha=0.3)
add_subplot_label(ax, '(a)')

# Panel (b): N content distribution
ax = axes[1]
ax.hist(orbital['HOMO_N_percent'], bins=10, alpha=0.7, label='HOMO', 
        color='#2ecc71', edgecolor='black')
ax.hist(orbital['LUMO_N_percent'], bins=10, alpha=0.7, label='LUMO',
        color='#e74c3c', edgecolor='black')
ax.set_xlabel('N Content (%)', fontsize=11)
ax.set_ylabel('Frequency', fontsize=11)
ax.set_title(f'N Content Distribution (n={len(orbital)})', 
             fontsize=11, fontweight='bold')
ax.legend(frameon=True)
ax.grid(axis='y', alpha=0.3)
add_subplot_label(ax, '(b)')

# Panel (c): Stacking energy distribution
ax = axes[2]
ax.hist(stack['interaction_energy_kcal_mol'], bins=8, color='#9b59b6', 
        alpha=0.7, edgecolor='black')
ax.axvline(x=stack['interaction_energy_kcal_mol'].mean(), color='red', 
           linestyle='--', linewidth=2, label='Mean')
ax.axvline(x=0, color='black', linestyle='-', linewidth=1, alpha=0.5)
ax.set_xlabel('Interaction Energy (kcal/mol)', fontsize=11)
ax.set_ylabel('Frequency', fontsize=11)
ax.set_title(f'Stacking Energy Distribution (n={len(stack)})', 
             fontsize=11, fontweight='bold')
ax.legend(frameon=True)
ax.grid(axis='y', alpha=0.3)
add_subplot_label(ax, '(c)')

plt.tight_layout()
save_figure(fig, 'statistical_summary_expanded')
plt.close()

print("\n" + "="*80)
print("FIGURE GENERATION COMPLETE")
print("="*80)
print(f"\n✓ 5 figures generated")
print(f"✓ 15 files created (5 figures × 3 formats)")
print(f"✓ All with scienceplots styling")
print(f"✓ All multi-panel figures labeled (a), (b), (c)")
print(f"✓ All saved at 300 dpi (PNG)")
print(f"\nOutput directory: {FIG_DIR}")
print("\nGenerated figures:")
print("  1. reorganization_energies_expanded.{png,pdf,pgf}")
print("  2. orbital_composition_expanded.{png,pdf,pgf}")
print("  3. stacking_energy_curves_expanded.{png,pdf,pgf}")
print("  4. property_correlations_expanded.{png,pdf,pgf}")
print("  5. statistical_summary_expanded.{png,pdf,pgf}")
