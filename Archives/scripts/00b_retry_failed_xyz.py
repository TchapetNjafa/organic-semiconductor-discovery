#!/usr/bin/env python3
"""
Retry failed XYZ generation with alternative methods
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path

def retry_with_alternatives(mol_id, smiles, output_dir):
    """Try multiple methods to generate geometry"""
    
    print(f"\nRetrying mol_id {mol_id}: {smiles}")
    
    # Method 1: Try with different random seed
    for seed in [42, 123, 456, 789, 2024]:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            mol = Chem.AddHs(mol)
            result = AllChem.EmbedMolecule(mol, randomSeed=seed)
            if result == 0:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
                save_xyz(mol, mol_id, smiles, output_dir)
                print(f"  ✓ Success with seed {seed}")
                return True
        except:
            continue
    
    # Method 2: Try with ETKDG (better for complex molecules)
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            mol = Chem.AddHs(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            result = AllChem.EmbedMolecule(mol, params)
            if result == 0:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
                save_xyz(mol, mol_id, smiles, output_dir)
                print(f"  ✓ Success with ETKDG")
                return True
    except Exception as e:
        print(f"  ETKDG failed: {e}")
    
    # Method 3: Try without optimization
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMultipleConfs(mol, numConfs=10, randomSeed=42)
            if mol.GetNumConformers() > 0:
                save_xyz(mol, mol_id, smiles, output_dir, conf_id=0)
                print(f"  ✓ Success with multiple conformers (no opt)")
                return True
    except Exception as e:
        print(f"  Multiple conformers failed: {e}")
    
    print(f"  ✗ All methods failed")
    return False

def save_xyz(mol, mol_id, smiles, output_dir, conf_id=-1):
    """Save molecule to XYZ file"""
    xyz_file = Path(output_dir) / f"{mol_id}.xyz"
    conf = mol.GetConformer(conf_id)
    
    with open(xyz_file, 'w') as f:
        f.write(f"{mol.GetNumAtoms()}\n")
        f.write(f"Generated from SMILES: {smiles}\n")
        
        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            symbol = atom.GetSymbol()
            f.write(f"{symbol:2s} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}\n")

if __name__ == '__main__':
    # Read original CSV
    pce_df = pd.read_csv('/home/tchapet/Post-Doc/ARTICLES DOSSIERS DES ARTICLES EN REDACTION/NOUVELS AXES DE RECHERCHE A REGARDER URGEMMENT/ARTICLES EN REDACTIONS/ARTICLE_MVOTO/MVOTOOPV/PUBCHEMQC_MULTIFONCTIONAL_2025/DATASET/dataset_pubchemqc_opv_17458.csv')
    
    # Read failed molecules
    failed_df = pd.read_csv('../data/failed_xyz_generation.csv')
    
    output_dir = Path('../data/input_geometries')
    
    print("="*60)
    print(f"Retrying {len(failed_df)} failed molecules")
    print("="*60)
    
    success_count = 0
    still_failed = []
    
    for _, row in failed_df.iterrows():
        mol_id = row['mol_id']
        
        # Get SMILES from original dataset
        mol_data = pce_df[pce_df['mol_id'] == mol_id]
        if len(mol_data) == 0:
            print(f"Mol {mol_id} not found in dataset")
            still_failed.append({'mol_id': mol_id, 'status': 'not_in_dataset'})
            continue
        
        smiles = mol_data.iloc[0]['SMILES']
        
        if retry_with_alternatives(mol_id, smiles, output_dir):
            success_count += 1
        else:
            still_failed.append({'mol_id': mol_id, 'status': 'all_methods_failed', 'smiles': smiles})
    
    print("\n" + "="*60)
    print("Retry Complete")
    print("="*60)
    print(f"Success: {success_count}/{len(failed_df)}")
    print(f"Still failed: {len(still_failed)}")
    
    if still_failed:
        still_failed_df = pd.DataFrame(still_failed)
        still_failed_df.to_csv('../data/still_failed_xyz.csv', index=False)
        print(f"\nStill failed molecules saved to: still_failed_xyz.csv")
        print("\nThese molecules may have:")
        print("  - Invalid SMILES")
        print("  - Highly strained structures")
        print("  - Unusual bonding patterns")
        print("\nYou can proceed with the remaining 17,446 molecules.")
