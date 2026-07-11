#!/usr/bin/env python3
"""
Utility functions for NTO analysis workflow
"""

import pandas as pd
import numpy as np
from pathlib import Path
import shutil

def load_pubchemqc_data():
    """Load the main PubChemQC dataset with PCE calculations"""
    csv_path = Path(__file__).parent.parent.parent / 'MVOTOOPV/PUBCHEMQC_MULTIFONCTIONAL_2025/DATASET/dataset_pubchemqc_opv_17458.csv'
    
    if not csv_path.exists():
        raise FileNotFoundError(f"Dataset not found: {csv_path}")
    
    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} molecules from PubChemQC dataset")
    return df

def prepare_xyz_geometries(df, output_dir, xyz_source_dir=None):
    """
    Prepare XYZ geometry files for NTO calculations
    
    Parameters:
    -----------
    df : DataFrame
        Molecules to process (should have 'mol_id' or 'gdb_id' column)
    output_dir : str/Path
        Directory to save XYZ files
    xyz_source_dir : str/Path
        Source directory containing XYZ files (if available)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if xyz_source_dir:
        # Copy existing XYZ files
        xyz_source = Path(xyz_source_dir)
        xyz_files = list(xyz_source.glob('*.xyz'))
        
        mol_ids = set(df['mol_id'].astype(str) if 'mol_id' in df else df['gdb_id'].astype(str))
        
        copied = 0
        for xyz_file in xyz_files:
            # Extract molecule ID from filename
            mol_id = xyz_file.stem.split('_')[0].replace('gdb ', '').replace('gdb', '')
            
            if mol_id in mol_ids:
                shutil.copy(xyz_file, output_dir / f"{mol_id}.xyz")
                copied += 1
        
        print(f"Copied {copied} XYZ files to {output_dir}")
        return copied
    
    else:
        print("No XYZ source directory provided. You'll need to generate geometries.")
        return 0

def filter_candidates_for_nto(df, criteria='opv_optimal'):
    """
    Filter molecules for NTO analysis based on criteria
    
    Parameters:
    -----------
    df : DataFrame
        Full dataset
    criteria : str
        'opv_optimal': PCE > 5%, SAScore < 5
        'top_pce': Top 1000 by PCE
        'top_combined': Top 1000 by PCE_SAScore
        'all': All molecules
    """
    if criteria == 'opv_optimal':
        filtered = df[(df['PCE_PCDTBT'] > 5) & (df['SAScore'] < 5)]
    elif criteria == 'top_pce':
        filtered = df.nlargest(1000, 'PCE_PCDTBT')
    elif criteria == 'top_combined':
        if 'PCE_SAScore' in df.columns:
            filtered = df.nlargest(1000, 'PCE_SAScore')
        else:
            print("PCE_SAScore not found, using top PCE instead")
            filtered = df.nlargest(1000, 'PCE_PCDTBT')
    else:  # 'all'
        filtered = df
    
    print(f"Filtered to {len(filtered)} molecules using '{criteria}' criteria")
    return filtered

def merge_nto_with_pce_data(nto_csv, pce_csv, output_csv):
    """Merge NTO results with existing PCE/SAScore data"""
    nto_df = pd.read_csv(nto_csv)
    pce_df = pd.read_csv(pce_csv)
    
    # Merge on molecule ID
    merged = pd.merge(
        pce_df, nto_df,
        left_on='mol_id' if 'mol_id' in pce_df else 'gdb_id',
        right_on='mol_id',
        how='left',
        suffixes=('', '_nto')
    )
    
    merged.to_csv(output_csv, index=False)
    print(f"Merged dataset saved to {output_csv}")
    print(f"Total molecules: {len(merged)}, with NTO data: {merged['status'].notna().sum()}")
    
    return merged

def calculate_nto_derived_metrics(df):
    """Calculate derived metrics from NTO properties"""
    
    # Exciton binding energy approximation
    if 'E_gap_fundamental' in df and 'S1_eV' in df:
        df['E_binding_eV'] = df['E_gap_fundamental'] - df['S1_eV']
    
    # Charge-transfer character classification
    if 'Sr_index' in df:
        df['CT_character'] = pd.cut(
            df['Sr_index'],
            bins=[0, 0.3, 0.6, 1.0],
            labels=['Strong CT', 'Moderate CT', 'Local']
        )
    
    # Single-excitation purity
    if 'Lambda_1' in df:
        df['single_exc_purity'] = df['Lambda_1'] > 0.8
    
    # OPV suitability score (0-1)
    score = 0
    if 'S1_eV' in df:
        # Optimal S1 range: 1.5-3.0 eV
        score += ((df['S1_eV'] >= 1.5) & (df['S1_eV'] <= 3.0)).astype(float) * 0.3
    if 'S1_f' in df:
        # Strong absorption: f > 0.1
        score += (df['S1_f'] > 0.1).astype(float) * 0.2
    if 'Sr_index' in df:
        # Moderate CT: 0.3 < Sr < 0.6
        score += ((df['Sr_index'] > 0.3) & (df['Sr_index'] < 0.6)).astype(float) * 0.3
    if 'D_index_Angstrom' in df:
        # Good separation: D > 2 Å
        score += (df['D_index_Angstrom'] > 2.0).astype(float) * 0.2
    
    df['NTO_OPV_score'] = score
    
    return df

def generate_summary_report(df, output_file):
    """Generate summary statistics report"""
    
    with open(output_file, 'w') as f:
        f.write("="*70 + "\n")
        f.write("NTO Analysis Summary Report\n")
        f.write("="*70 + "\n\n")
        
        f.write(f"Total molecules analyzed: {len(df)}\n")
        f.write(f"Successful calculations: {(df['status'] == 'success').sum()}\n\n")
        
        if 'S1_eV' in df.columns:
            f.write("Excitation Energy (S1):\n")
            f.write(f"  Mean: {df['S1_eV'].mean():.2f} ± {df['S1_eV'].std():.2f} eV\n")
            f.write(f"  Range: {df['S1_eV'].min():.2f} - {df['S1_eV'].max():.2f} eV\n")
            f.write(f"  OPV optimal (1.5-3.0 eV): {((df['S1_eV'] >= 1.5) & (df['S1_eV'] <= 3.0)).sum()}\n\n")
        
        if 'S1_f' in df.columns:
            f.write("Oscillator Strength (S1):\n")
            f.write(f"  Mean: {df['S1_f'].mean():.3f} ± {df['S1_f'].std():.3f}\n")
            f.write(f"  Strong absorption (f > 0.1): {(df['S1_f'] > 0.1).sum()}\n\n")
        
        if 'Sr_index' in df.columns:
            f.write("Electron-Hole Overlap (Sr):\n")
            f.write(f"  Mean: {df['Sr_index'].mean():.3f} ± {df['Sr_index'].std():.3f}\n")
            f.write(f"  Moderate CT (0.3-0.6): {((df['Sr_index'] > 0.3) & (df['Sr_index'] < 0.6)).sum()}\n\n")
        
        if 'D_index_Angstrom' in df.columns:
            f.write("Charge-Transfer Distance (D):\n")
            f.write(f"  Mean: {df['D_index_Angstrom'].mean():.2f} ± {df['D_index_Angstrom'].std():.2f} Å\n")
            f.write(f"  Good separation (> 2 Å): {(df['D_index_Angstrom'] > 2.0).sum()}\n\n")
        
        if 'NTO_OPV_score' in df.columns:
            f.write("NTO-based OPV Suitability Score:\n")
            f.write(f"  Mean: {df['NTO_OPV_score'].mean():.2f} ± {df['NTO_OPV_score'].std():.2f}\n")
            f.write(f"  High score (> 0.7): {(df['NTO_OPV_score'] > 0.7).sum()}\n\n")
    
    print(f"Summary report saved to {output_file}")

if __name__ == '__main__':
    # Test utilities
    print("Testing utility functions...")
    
    try:
        df = load_pubchemqc_data()
        print(f"✓ Loaded {len(df)} molecules")
        print(f"  Columns: {', '.join(df.columns[:10])}...")
    except Exception as e:
        print(f"✗ Error loading data: {e}")
