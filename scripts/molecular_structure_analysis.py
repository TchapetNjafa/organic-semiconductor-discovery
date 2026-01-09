#!/usr/bin/env python3
"""
Molecular Structure Analysis and Visualization
===============================================
This script analyzes the chemical structures of top candidates and creates
annotated molecular structure figures with functional group identification.

Author: Analysis for Computational Materials Science submission
Date: 2026-01-08
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw, Fragments
from rdkit.Chem.Draw import IPythonConsole, rdMolDraw2D
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300

# Top 7 molecule IDs from the manuscript
TOP_MOLECULES = {
    977: 'Acceptor',
    1712: 'Acceptor',
    4550: 'Donor',
    7801: 'Acceptor',
    11029: 'Donor',
    17851: 'Donor',
    20778: 'Donor'
}

def load_molecule_data():
    """Load SMILES and properties for top molecules"""
    # Load from predictions file
    df_pred = pd.read_csv('DATASET/predictions_molecules_cibles.csv')
    
    # Extract molecule IDs and SMILES
    molecules = []
    for idx, row in df_pred.iterrows():
        mol_id_str = row['mol_id']
        # Extract numeric ID from string like "Accepteur_977" or "Donneur_4550"
        mol_id = int(mol_id_str.split('_')[1])
        
        molecules.append({
            'mol_id': mol_id,
            'type': TOP_MOLECULES.get(mol_id, 'Unknown'),
            'SMILES': row['SMILES'],
            'MW': row['MW'],
            'LogP': row['LogP'],
            'HBA': row['HBA'],
            'HBD': row['HBD'],
            'TPSA': row['TPSA'],
            'QED': row['QED Weighted']
        })
    
    return pd.DataFrame(molecules)

def identify_functional_groups(mol):
    """Identify key functional groups in a molecule"""
    if mol is None:
        return []
    
    groups = []
    
    # Aromatic rings
    if Fragments.fr_benzene(mol) > 0:
        groups.append(f"Benzene rings ({Fragments.fr_benzene(mol)})")
    
    # Heteroaromatic
    if Fragments.fr_thiophene(mol) > 0:
        groups.append(f"Thiophene ({Fragments.fr_thiophene(mol)})")
    if Fragments.fr_furan(mol) > 0:
        groups.append(f"Furan ({Fragments.fr_furan(mol)})")
    if Fragments.fr_pyridine(mol) > 0:
        groups.append(f"Pyridine ({Fragments.fr_pyridine(mol)})")
    
    # Functional groups
    if Fragments.fr_C_O(mol) > 0:
        groups.append(f"C=O ({Fragments.fr_C_O(mol)})")
    if Fragments.fr_COO(mol) > 0:
        groups.append(f"Carboxyl ({Fragments.fr_COO(mol)})")
    if Fragments.fr_NH2(mol) > 0:
        groups.append(f"Amino ({Fragments.fr_NH2(mol)})")
    if Fragments.fr_Ar_OH(mol) > 0:
        groups.append(f"Phenol ({Fragments.fr_Ar_OH(mol)})")
    if Fragments.fr_Ar_N(mol) > 0:
        groups.append(f"Aromatic N ({Fragments.fr_Ar_N(mol)})")
    if Fragments.fr_nitrile(mol) > 0:
        groups.append(f"Nitrile ({Fragments.fr_nitrile(mol)})")
    if Fragments.fr_azide(mol) > 0:
        groups.append(f"Azide ({Fragments.fr_azide(mol)})")
    
    # Count heteroatoms
    n_atoms = mol.GetNumAtoms()
    heteroatoms = [atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetSymbol() not in ['C', 'H']]
    if heteroatoms:
        hetero_counts = Counter(heteroatoms)
        groups.append(f"Heteroatoms: {dict(hetero_counts)}")
    
    return groups

def analyze_conjugation(mol):
    """Analyze conjugation patterns"""
    if mol is None:
        return {}
    
    # Count aromatic atoms
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    total_atoms = mol.GetNumHeavyAtoms()
    
    # Count aromatic bonds
    aromatic_bonds = sum(1 for bond in mol.GetBonds() if bond.GetIsAromatic())
    total_bonds = mol.GetNumBonds()
    
    # Count rotatable bonds (measure of rigidity/planarity)
    rotatable = Descriptors.NumRotatableBonds(mol)
    
    return {
        'aromatic_atoms': aromatic_atoms,
        'aromatic_fraction': aromatic_atoms / total_atoms if total_atoms > 0 else 0,
        'aromatic_bonds': aromatic_bonds,
        'rotatable_bonds': rotatable,
        'rigidity_score': 1 - (rotatable / total_bonds) if total_bonds > 0 else 0
    }

def create_structure_table(df_mols):
    """Create a comprehensive structure-property table"""
    results = []
    
    for idx, row in df_mols.iterrows():
        mol = Chem.MolFromSmiles(row['SMILES'])
        
        if mol is not None:
            # Get functional groups
            func_groups = identify_functional_groups(mol)
            
            # Get conjugation info
            conj_info = analyze_conjugation(mol)
            
            results.append({
                'Mol_ID': row['mol_id'],
                'Type': row['type'],
                'Formula': Chem.rdMolDescriptors.CalcMolFormula(mol),
                'Heavy_Atoms': mol.GetNumHeavyAtoms(),
                'Aromatic_Fraction': f"{conj_info['aromatic_fraction']:.2f}",
                'Rotatable_Bonds': conj_info['rotatable_bonds'],
                'Functional_Groups': '; '.join(func_groups[:3]),  # Top 3
                'LogP': f"{row['LogP']:.2f}",
                'TPSA': f"{row['TPSA']:.1f}",
                'QED': f"{row['QED']:.3f}"
            })
    
    return pd.DataFrame(results)

def draw_molecules_grid(df_mols, output_file='figures/top_molecules_structures.png'):
    """Draw all molecules in a grid with annotations"""
    mols = []
    legends = []
    
    for idx, row in df_mols.iterrows():
        mol = Chem.MolFromSmiles(row['SMILES'])
        if mol is not None:
            AllChem.Compute2DCoords(mol)
            mols.append(mol)
            legends.append(f"Mol {row['mol_id']} ({row['type']})")
    
    # Create grid image
    img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(400, 400),
                                legends=legends, returnPNG=False)
    
    img.save(output_file, dpi=(300, 300))
    print(f"✓ Saved molecular structures to: {output_file}")
    
    return img

def main():
    """Main execution"""
    print("="*60)
    print("Molecular Structure Analysis")
    print("="*60)
    
    # Load data
    print("\n📊 Loading molecule data...")
    df_mols = load_molecule_data()
    print(f"✓ Loaded {len(df_mols)} molecules")
    print(df_mols[['mol_id', 'type', 'SMILES']])
    
    # Create structure-property table
    print("\n🔬 Analyzing chemical structures...")
    df_analysis = create_structure_table(df_mols)
    
    # Save table
    df_analysis.to_csv('DATASET/molecular_structure_analysis.csv', index=False)
    print(f"✓ Saved analysis to: DATASET/molecular_structure_analysis.csv")
    
    # Print table
    print("\n" + "="*60)
    print("STRUCTURE-PROPERTY ANALYSIS")
    print("="*60)
    print(df_analysis.to_string(index=False))
    
    # Draw molecules
    print("\n🎨 Creating molecular structure visualizations...")
    draw_molecules_grid(df_mols)
    
    print("\n✅ Analysis complete!")
    print("="*60)

if __name__ == "__main__":
    main()

