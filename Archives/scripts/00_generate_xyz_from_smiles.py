#!/usr/bin/env python3
"""
Generate XYZ files from SMILES in dataset_pubchemqc_opv_17458.csv
Creates 3D geometries for all 17,458 molecules
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path
from joblib import Parallel, delayed
import multiprocessing

def smiles_to_xyz(mol_id, smiles, output_dir):
    """Convert SMILES to XYZ file with 3D coordinates"""
    try:
        # Create molecule from SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'mol_id': mol_id, 'status': 'invalid_smiles'}
        
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates
        result = AllChem.EmbedMolecule(mol, randomSeed=42)
        if result != 0:
            return {'mol_id': mol_id, 'status': 'embed_failed'}
        
        # Optimize geometry with UFF
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        
        # Write XYZ file
        xyz_file = Path(output_dir) / f"{mol_id}.xyz"
        conf = mol.GetConformer()
        
        with open(xyz_file, 'w') as f:
            f.write(f"{mol.GetNumAtoms()}\n")
            f.write(f"Generated from SMILES: {smiles}\n")
            
            for i, atom in enumerate(mol.GetAtoms()):
                pos = conf.GetAtomPosition(i)
                symbol = atom.GetSymbol()
                f.write(f"{symbol:2s} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}\n")
        
        return {'mol_id': mol_id, 'status': 'success'}
        
    except Exception as e:
        return {'mol_id': mol_id, 'status': f'error: {str(e)}'}

def generate_all_xyz(csv_file, output_dir, n_jobs=-1):
    """Generate XYZ files for all molecules in parallel"""
    
    # Read CSV
    df = pd.read_csv(csv_file)
    print(f"Found {len(df)} molecules in {csv_file}")
    
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Check for existing files
    existing = set(f.stem for f in output_dir.glob('*.xyz'))
    to_process = df[~df['mol_id'].astype(str).isin(existing)]
    
    print(f"Already have {len(existing)} XYZ files")
    print(f"Need to generate {len(to_process)} XYZ files")
    
    if len(to_process) == 0:
        print("All XYZ files already exist!")
        return
    
    # Parallel processing
    print(f"Using {n_jobs if n_jobs > 0 else multiprocessing.cpu_count()} CPU cores")
    
    results = Parallel(n_jobs=n_jobs, verbose=10)(
        delayed(smiles_to_xyz)(row['mol_id'], row['SMILES'], output_dir)
        for _, row in to_process.iterrows()
    )
    
    # Summary
    results_df = pd.DataFrame(results)
    success = (results_df['status'] == 'success').sum()
    failed = len(results_df) - success
    
    print(f"\n{'='*60}")
    print(f"XYZ Generation Complete")
    print(f"{'='*60}")
    print(f"Success: {success}")
    print(f"Failed: {failed}")
    print(f"Total XYZ files: {len(list(output_dir.glob('*.xyz')))}")
    
    # Save failed molecules
    if failed > 0:
        failed_df = results_df[results_df['status'] != 'success']
        failed_df.to_csv(output_dir.parent / 'failed_xyz_generation.csv', index=False)
        print(f"Failed molecules saved to: failed_xyz_generation.csv")

if __name__ == '__main__':
    CSV_FILE = '/home/tchapet/Post-Doc/ARTICLES DOSSIERS DES ARTICLES EN REDACTION/NOUVELS AXES DE RECHERCHE A REGARDER URGEMMENT/ARTICLES EN REDACTIONS/ARTICLE_MVOTO/MVOTOOPV/PUBCHEMQC_MULTIFONCTIONAL_2025/DATASET/dataset_pubchemqc_opv_17458.csv'
    OUTPUT_DIR = '../data/input_geometries'
    
    print("="*60)
    print("Generating XYZ files from SMILES")
    print("="*60)
    print(f"Input: {CSV_FILE}")
    print(f"Output: {OUTPUT_DIR}")
    print("")
    
    generate_all_xyz(CSV_FILE, OUTPUT_DIR, n_jobs=-1)
