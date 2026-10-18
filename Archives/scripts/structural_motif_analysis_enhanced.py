#!/usr/bin/env python3
"""
Enhanced Structural Motif Analysis Script
==========================================

Analyzes structural patterns for top N molecules (100, 500, or 1000)
to identify donor-acceptor patterns and structure-property relationships.

Usage: python structural_motif_analysis_enhanced.py [N_MOLECULES]
Example: python structural_motif_analysis_enhanced.py 100

Author: Enhanced for ARTICLE_MVOTO project
Date: 2026-01-09
Phase: 4.2 - Extended Structural Analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Fragments, Lipinski
from rdkit.Chem import Draw
from collections import Counter
import warnings
import sys
import time
warnings.filterwarnings('ignore')

# Get number of molecules from command line
N_MOLECULES = int(sys.argv[1]) if len(sys.argv) > 1 else 100

# Set style
plt.style.use('seaborn-v0_8-paper')
sns.set_palette("husl")

# Define paths
DATASET_DIR = Path("DATASET")
FIGURES_DIR = Path("figures")
FIGURES_DIR.mkdir(exist_ok=True)

print("="*80)
print(f"ENHANCED STRUCTURAL MOTIF ANALYSIS - TOP {N_MOLECULES} MOLECULES")
print("="*80)
print(f"Start time: {time.strftime('%Y-%m-%d %H:%M:%S')}")

# ============================================================================
# SECTION 1: Load Data
# ============================================================================
print("\nSECTION 1: Loading Data")
print("-"*80)

df_pce = pd.read_csv(DATASET_DIR / "dataset_pubchemqc_opv_17458.csv")
df_candidates = pd.read_csv(DATASET_DIR / "predictions_molecules_cibles.csv")

print(f"✓ Loaded PCE data: {len(df_pce)} molecules")
print(f"✓ Loaded top candidates: {len(df_candidates)} molecules")

# ============================================================================
# SECTION 2: Functional Group Analysis Function
# ============================================================================

def analyze_functional_groups(smiles):
    """Comprehensive functional group analysis"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    groups = {
        # Aromatic systems
        'Benzene_rings': Fragments.fr_benzene(mol),
        'Aromatic_rings': Lipinski.NumAromaticRings(mol),
        'Heteroaromatic_rings': Lipinski.NumAromaticHeterocycles(mol),
        
        # Nitrogen-containing groups
        'Amines': Fragments.fr_NH2(mol) + Fragments.fr_NH1(mol) + Fragments.fr_NH0(mol),
        'Nitriles': Fragments.fr_nitrile(mol),
        'Nitro_groups': Fragments.fr_nitro(mol),
        'Azides': Fragments.fr_azide(mol),
        'Imines': Fragments.fr_Imine(mol),
        
        # Oxygen-containing groups
        'Carbonyls': Fragments.fr_C_O(mol),
        'Esters': Fragments.fr_ester(mol),
        'Ethers': Fragments.fr_ether(mol),
        'Alcohols': Fragments.fr_Al_OH(mol),
        'Phenols': Fragments.fr_phenol(mol),
        
        # Sulfur-containing groups
        'Thiophenes': Fragments.fr_thiophene(mol),
        'Sulfones': Fragments.fr_sulfone(mol),
        'Sulfonamides': Fragments.fr_sulfonamd(mol),
        
        # Conjugation features
        'Conjugated_double_bonds': len([bond for bond in mol.GetBonds() 
                                       if bond.GetBondType() == Chem.BondType.DOUBLE 
                                       and bond.GetIsConjugated()]),
        
        # Heteroatoms
        'N_atoms': sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'N'),
        'O_atoms': sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O'),
        'S_atoms': sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'S'),
        'F_atoms': sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'F'),
        'Cl_atoms': sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'Cl'),
        
        # Structural features
        'Rotatable_bonds': Lipinski.NumRotatableBonds(mol),
        'HBD': Lipinski.NumHDonors(mol),
        'HBA': Lipinski.NumHAcceptors(mol),
        'Aromatic_fraction': Descriptors.NumAromaticRings(mol) / max(1, Descriptors.RingCount(mol)) if Descriptors.RingCount(mol) > 0 else 0,
    }
    
    return groups

# ============================================================================
# SECTION 3: Analyze Top N Molecules
# ============================================================================
print(f"\nSECTION 2: Analyzing Top {N_MOLECULES} Molecules by PCE_SAScore")
print("-"*80)

# Calculate PCE_SAScore and sort
df_pce['PCE_SAScore_PCDTBT'] = df_pce['pce_pcdtbt(%)'] - df_pce['sas1(%)']
df_top = df_pce.nlargest(N_MOLECULES, 'PCE_SAScore_PCDTBT').copy()

print(f"✓ Selected top {len(df_top)} molecules")
print(f"  PCE_SAScore range: {df_top['PCE_SAScore_PCDTBT'].min():.2f} to {df_top['PCE_SAScore_PCDTBT'].max():.2f}")

# Analyze functional groups
print(f"\nAnalyzing functional groups (this may take a few minutes)...")
start_time = time.time()

df_analysis_full = []
failed = 0

for idx, row in df_top.iterrows():
    if (idx + 1) % 50 == 0:
        elapsed = time.time() - start_time
        rate = (idx + 1) / elapsed
        remaining = (N_MOLECULES - idx - 1) / rate
        print(f"  Progress: {idx+1}/{N_MOLECULES} ({100*(idx+1)/N_MOLECULES:.1f}%) - ETA: {remaining:.0f}s")
    
    smiles = row['SMILES']
    groups = analyze_functional_groups(smiles)
    
    if groups is None:
        failed += 1
        continue
    
    groups.update({
        'mol_id': row['mol_id'],
        'SMILES': smiles,
        'HOMO': row['HOMO(eV)'],
        'LUMO': row['LUMO(eV)'],
        'GAP': row['GAP(eV)'],
        'PCE_PCBM': row['pce_pcbm(%)'],
        'PCE_PCDTBT': row['pce_pcdtbt(%)'],
        'SAScore': row['sas1(%)'],
        'PCE_SAScore': row['PCE_SAScore_PCDTBT']
    })
    df_analysis_full.append(groups)

df_full = pd.DataFrame(df_analysis_full)

elapsed = time.time() - start_time
print(f"\n✓ Analysis complete in {elapsed:.1f}s")
print(f"  Successful: {len(df_full)}/{N_MOLECULES}")
print(f"  Failed: {failed}")

# Save results
output_file = DATASET_DIR / f"structural_motif_analysis_top{N_MOLECULES}.csv"
df_full.to_csv(output_file, index=False)
print(f"✓ Saved results to: {output_file}")

# ============================================================================
# SECTION 4: Statistical Analysis
# ============================================================================
print("\n" + "="*80)
print("SECTION 3: Statistical Analysis")
print("="*80)

# Correlations with PCE
print("\nCorrelations with PCE (PCDTBT):")
print("-"*80)

structural_features = ['Aromatic_rings', 'Heteroaromatic_rings', 'N_atoms', 'O_atoms', 
                      'S_atoms', 'Conjugated_double_bonds', 'Rotatable_bonds', 
                      'Aromatic_fraction', 'HBD', 'HBA']

correlations = []
for feature in structural_features:
    corr = df_full[[feature, 'PCE_PCDTBT']].corr().iloc[0, 1]
    correlations.append({'Feature': feature, 'Correlation': corr})
    print(f"  {feature:30s}: {corr:+.3f}")

df_correlations = pd.DataFrame(correlations)

# Summary statistics
print("\nSummary Statistics:")
print("-"*80)
summary_features = ['Aromatic_rings', 'N_atoms', 'O_atoms', 'S_atoms', 
                   'Conjugated_double_bonds', 'Rotatable_bonds']

for feature in summary_features:
    print(f"  {feature:30s}: mean={df_full[feature].mean():.2f}, std={df_full[feature].std():.2f}, max={df_full[feature].max():.0f}")

# ============================================================================
# SECTION 5: Identify Donor/Acceptor Patterns
# ============================================================================
print("\n" + "="*80)
print("SECTION 4: Donor/Acceptor Pattern Classification")
print("="*80)

# Classify based on HOMO-LUMO
# Donors: higher HOMO (electron-rich)
# Acceptors: lower LUMO (electron-deficient)

median_homo = df_full['HOMO'].median()
median_lumo = df_full['LUMO'].median()

df_full['Predicted_Type'] = 'Intermediate'
df_full.loc[df_full['HOMO'] > median_homo, 'Predicted_Type'] = 'Donor-like'
df_full.loc[df_full['LUMO'] < median_lumo, 'Predicted_Type'] = 'Acceptor-like'

print(f"\nClassification based on HOMO/LUMO:")
print(f"  Donor-like (HOMO > {median_homo:.2f} eV): {(df_full['Predicted_Type']=='Donor-like').sum()}")
print(f"  Acceptor-like (LUMO < {median_lumo:.2f} eV): {(df_full['Predicted_Type']=='Acceptor-like').sum()}")
print(f"  Intermediate: {(df_full['Predicted_Type']=='Intermediate').sum()}")

# Compare structural features
donors = df_full[df_full['Predicted_Type'] == 'Donor-like']
acceptors = df_full[df_full['Predicted_Type'] == 'Acceptor-like']

print(f"\nDonor-like molecules (n={len(donors)}):")
print("-"*40)
print(f"  Avg N atoms: {donors['N_atoms'].mean():.2f}")
print(f"  Avg O atoms: {donors['O_atoms'].mean():.2f}")
print(f"  Avg aromatic rings: {donors['Aromatic_rings'].mean():.2f}")
print(f"  Avg conjugated bonds: {donors['Conjugated_double_bonds'].mean():.2f}")

print(f"\nAcceptor-like molecules (n={len(acceptors)}):")
print("-"*40)
print(f"  Avg N atoms: {acceptors['N_atoms'].mean():.2f}")
print(f"  Avg O atoms: {acceptors['O_atoms'].mean():.2f}")
print(f"  Avg aromatic rings: {acceptors['Aromatic_rings'].mean():.2f}")
print(f"  Avg conjugated bonds: {acceptors['Conjugated_double_bonds'].mean():.2f}")

# ============================================================================
# SECTION 6: Visualization
# ============================================================================
print("\n" + "="*80)
print("SECTION 5: Creating Visualizations")
print("="*80)

# Figure 1: Correlation heatmap
fig1, ax1 = plt.subplots(figsize=(10, 8))
corr_matrix = df_full[structural_features + ['PCE_PCDTBT', 'SAScore']].corr()
sns.heatmap(corr_matrix, annot=True, fmt='.2f', cmap='coolwarm', center=0, 
            square=True, linewidths=1, cbar_kws={"shrink": 0.8}, ax=ax1)
ax1.set_title(f'Structure-Property Correlations (Top {N_MOLECULES} Molecules)', 
              fontsize=14, fontweight='bold', pad=20)
plt.tight_layout()
fig1.savefig(FIGURES_DIR / f'correlation_heatmap_top{N_MOLECULES}.pdf', dpi=300, bbox_inches='tight')
fig1.savefig(FIGURES_DIR / f'correlation_heatmap_top{N_MOLECULES}.png', dpi=300, bbox_inches='tight')
print(f"✓ Saved: correlation_heatmap_top{N_MOLECULES}.pdf/png")

# Figure 2: Donor vs Acceptor comparison
fig2, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: Heteroatom content
ax = axes[0, 0]
features = ['N_atoms', 'O_atoms', 'S_atoms']
x = np.arange(len(features))
width = 0.35
donor_vals = [donors[f].mean() for f in features]
acceptor_vals = [acceptors[f].mean() for f in features]
ax.bar(x - width/2, donor_vals, width, label='Donor-like', alpha=0.8, edgecolor='black')
ax.bar(x + width/2, acceptor_vals, width, label='Acceptor-like', alpha=0.8, edgecolor='black')
ax.set_ylabel('Average Count', fontweight='bold')
ax.set_title('Heteroatom Content', fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(features)
ax.legend()
ax.grid(axis='y', alpha=0.3)

# Plot 2: Aromatic features
ax = axes[0, 1]
features = ['Aromatic_rings', 'Heteroaromatic_rings', 'Benzene_rings']
x = np.arange(len(features))
donor_vals = [donors[f].mean() for f in features]
acceptor_vals = [acceptors[f].mean() for f in features]
ax.bar(x - width/2, donor_vals, width, label='Donor-like', alpha=0.8, edgecolor='black')
ax.bar(x + width/2, acceptor_vals, width, label='Acceptor-like', alpha=0.8, edgecolor='black')
ax.set_ylabel('Average Count', fontweight='bold')
ax.set_title('Aromatic Features', fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels([f.replace('_', '\n') for f in features], fontsize=9)
ax.legend()
ax.grid(axis='y', alpha=0.3)

# Plot 3: Functional groups
ax = axes[1, 0]
features = ['Nitriles', 'Azides', 'Carbonyls', 'Phenols']
x = np.arange(len(features))
donor_vals = [donors[f].mean() for f in features]
acceptor_vals = [acceptors[f].mean() for f in features]
ax.bar(x - width/2, donor_vals, width, label='Donor-like', alpha=0.8, edgecolor='black')
ax.bar(x + width/2, acceptor_vals, width, label='Acceptor-like', alpha=0.8, edgecolor='black')
ax.set_ylabel('Average Count', fontweight='bold')
ax.set_title('Key Functional Groups', fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(features)
ax.legend()
ax.grid(axis='y', alpha=0.3)

# Plot 4: Structural flexibility
ax = axes[1, 1]
features = ['Conjugated_double_bonds', 'Rotatable_bonds']
x = np.arange(len(features))
donor_vals = [donors[f].mean() for f in features]
acceptor_vals = [acceptors[f].mean() for f in features]
ax.bar(x - width/2, donor_vals, width, label='Donor-like', alpha=0.8, edgecolor='black')
ax.bar(x + width/2, acceptor_vals, width, label='Acceptor-like', alpha=0.8, edgecolor='black')
ax.set_ylabel('Average Count', fontweight='bold')
ax.set_title('Conjugation & Flexibility', fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels([f.replace('_', '\n') for f in features], fontsize=9)
ax.legend()
ax.grid(axis='y', alpha=0.3)

plt.suptitle(f'Donor-like vs Acceptor-like Structural Patterns (Top {N_MOLECULES} Molecules)', 
             fontsize=16, fontweight='bold', y=1.00)
plt.tight_layout()
fig2.savefig(FIGURES_DIR / f'donor_acceptor_patterns_top{N_MOLECULES}.pdf', dpi=300, bbox_inches='tight')
fig2.savefig(FIGURES_DIR / f'donor_acceptor_patterns_top{N_MOLECULES}.png', dpi=300, bbox_inches='tight')
print(f"✓ Saved: donor_acceptor_patterns_top{N_MOLECULES}.pdf/png")

# ============================================================================
# SECTION 7: Summary Report
# ============================================================================
print("\n" + "="*80)
print("SUMMARY REPORT")
print("="*80)

report = f"""
ENHANCED STRUCTURAL MOTIF ANALYSIS REPORT
==========================================

Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}
Molecules Analyzed: {N_MOLECULES}
Successful: {len(df_full)}
Failed: {failed}

KEY FINDINGS:
-------------

1. Donor-like molecules (n={len(donors)}):
   - Higher HOMO (> {median_homo:.2f} eV)
   - Avg N atoms: {donors['N_atoms'].mean():.2f}
   - Avg O atoms: {donors['O_atoms'].mean():.2f}
   - Avg aromatic rings: {donors['Aromatic_rings'].mean():.2f}

2. Acceptor-like molecules (n={len(acceptors)}):
   - Lower LUMO (< {median_lumo:.2f} eV)
   - Avg N atoms: {acceptors['N_atoms'].mean():.2f}
   - Avg O atoms: {acceptors['O_atoms'].mean():.2f}
   - Avg aromatic rings: {acceptors['Aromatic_rings'].mean():.2f}

3. Top correlations with PCE:
{df_correlations.nlargest(5, 'Correlation', keep='all').to_string(index=False)}

FILES GENERATED:
----------------
- {output_file}
- correlation_heatmap_top{N_MOLECULES}.pdf/png
- donor_acceptor_patterns_top{N_MOLECULES}.pdf/png

NEXT STEPS:
-----------
1. Review figures in {FIGURES_DIR}/
2. Analyze detailed results in {output_file}
3. Update manuscript with new findings
4. Consider running with larger N (500 or 1000) for more robust statistics
"""

print(report)

# Save report
report_file = DATASET_DIR / f"structural_motif_analysis_top{N_MOLECULES}_report.txt"
with open(report_file, 'w') as f:
    f.write(report)
print(f"\n✓ Saved report to: {report_file}")

print("\n" + "="*80)
print("ANALYSIS COMPLETE!")
print("="*80)
print(f"Total time: {time.time() - start_time:.1f}s")
