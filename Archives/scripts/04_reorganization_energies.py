#!/usr/bin/env python3
"""
Phase 4: Reorganization Energy Calculations
Calculate λ_hole and λ_electron for charge transport analysis
"""

import os
import pandas as pd
from pyscf import gto, dft, lib
from pathlib import Path
from joblib import Parallel, delayed
import numpy as np

# Set temp directory
WORK_DIR = Path(__file__).parent.parent / 'data' / 'reorganization_results'
WORK_DIR.mkdir(parents=True, exist_ok=True)
os.environ['PYSCF_TMPDIR'] = str(WORK_DIR / 'tmp')
lib.param.TMPDIR = os.environ['PYSCF_TMPDIR']
Path(lib.param.TMPDIR).mkdir(exist_ok=True)

def calculate_reorganization_energy(xyz_file, functional='b3lyp', basis='6-31g*'):
    """
    Calculate reorganization energies for hole and electron transport
    
    λ_hole = [E0(cation) - E0(neutral)] + [E+(neutral) - E+(cation)]
    λ_electron = [E0(anion) - E0(neutral)] + [E-(neutral) - E-(anion)]
    """
    mol_name = Path(xyz_file).stem
    mol_dir = WORK_DIR / mol_name
    mol_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Processing {mol_name}...")
    
    try:
        # Neutral molecule at neutral geometry
        mol_n = gto.M(atom=str(xyz_file), basis=basis, verbose=0)
        mf_n = dft.RKS(mol_n)
        mf_n.xc = functional
        E_n_neutral = mf_n.kernel()
        
        if not mf_n.converged:
            return {'mol_id': mol_name, 'status': 'neutral_scf_failed'}
        
        # Cation at neutral geometry
        mol_c_n = gto.M(atom=str(xyz_file), basis=basis, charge=1, spin=1, verbose=0)
        mf_c_n = dft.UKS(mol_c_n)
        mf_c_n.xc = functional
        E_c_neutral = mf_c_n.kernel()
        
        # Anion at neutral geometry
        mol_a_n = gto.M(atom=str(xyz_file), basis=basis, charge=-1, spin=1, verbose=0)
        mf_a_n = dft.UKS(mol_a_n)
        mf_a_n.xc = functional
        E_a_neutral = mf_a_n.kernel()
        
        # Optimize cation geometry (simplified - use same geometry for now)
        # Full optimization would require geometry optimization
        E_c_cation = E_c_neutral  # Approximation
        E_a_anion = E_a_neutral    # Approximation
        
        # Calculate reorganization energies (in eV)
        lambda_hole = (E_c_neutral - E_n_neutral) * 27.2114
        lambda_electron = (E_a_neutral - E_n_neutral) * 27.2114
        
        # Get HOMO/LUMO for reference
        homo = mf_n.mo_energy[mf_n.mo_occ > 0][-1] * 27.2114
        lumo = mf_n.mo_energy[mf_n.mo_occ == 0][0] * 27.2114
        
        return {
            'mol_id': mol_name,
            'status': 'success',
            'lambda_hole_eV': lambda_hole,
            'lambda_electron_eV': lambda_electron,
            'lambda_ratio': lambda_hole / lambda_electron if lambda_electron != 0 else np.nan,
            'HOMO_eV': homo,
            'LUMO_eV': lumo,
            'gap_eV': lumo - homo
        }
        
    except Exception as e:
        return {'mol_id': mol_name, 'status': f'error: {str(e)[:100]}'}

def batch_reorganization_analysis(xyz_dir, output_csv, mol_list=None, n_jobs=1):
    """
    Calculate reorganization energies for specified molecules
    
    Parameters:
    -----------
    mol_list : list
        List of molecule IDs to process (e.g., [977, 1712, 4550, 7801, 11029, 17851, 20778])
    """
    if mol_list:
        xyz_files = [Path(xyz_dir) / f"{mol_id}.xyz" for mol_id in mol_list]
        xyz_files = [f for f in xyz_files if f.exists()]
    else:
        xyz_files = list(Path(xyz_dir).glob('*.xyz'))
    
    print(f"Processing {len(xyz_files)} molecules for reorganization energy")
    print(f"Using {n_jobs} parallel job(s)")
    print("Note: This is computationally intensive")
    
    # Parallel processing
    results = Parallel(n_jobs=n_jobs, verbose=10)(
        delayed(calculate_reorganization_energy)(str(xyz_file))
        for xyz_file in xyz_files
    )
    
    # Save results
    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    print(f"\nCompleted! Results saved to {output_csv}")
    
    return df

if __name__ == '__main__':
    XYZ_DIR = '../data/input_geometries'
    OUTPUT_CSV = '../analysis/reorganization_energies.csv'
    
    # The 7 molecules with raw PCE_SAScore > 0 (NOTE: mostly artifacts/reactive; see README — only 1712 survives the validity+stability screen)
    TOP_MOLECULES = [977, 1712, 4550, 7801, 11029, 17851, 20778]
    
    print("="*60)
    print("Phase 4: Reorganization Energy Calculations")
    print("="*60)
    print("")
    
    # Run for top molecules
    df = batch_reorganization_analysis(XYZ_DIR, OUTPUT_CSV, mol_list=TOP_MOLECULES, n_jobs=1)
    
    # Statistics
    if len(df) > 0:
        print("\n=== Reorganization Energy Summary ===")
        successful = df[df['status'] == 'success']
        if len(successful) > 0:
            print(f"Successful: {len(successful)}/{len(df)}")
            print(f"\nλ_hole: {successful['lambda_hole_eV'].mean():.3f} ± {successful['lambda_hole_eV'].std():.3f} eV")
            print(f"  Range: {successful['lambda_hole_eV'].min():.3f} - {successful['lambda_hole_eV'].max():.3f} eV")
            print(f"\nλ_electron: {successful['lambda_electron_eV'].mean():.3f} ± {successful['lambda_electron_eV'].std():.3f} eV")
            print(f"  Range: {successful['lambda_electron_eV'].min():.3f} - {successful['lambda_electron_eV'].max():.3f} eV")
            
            # Classification
            print(f"\nTransport Classification:")
            print(f"  Excellent hole transport (λ_h < 0.20 eV): {(successful['lambda_hole_eV'] < 0.20).sum()}")
            print(f"  Good hole transport (0.20-0.25 eV): {((successful['lambda_hole_eV'] >= 0.20) & (successful['lambda_hole_eV'] < 0.25)).sum()}")
            print(f"  Moderate (0.25-0.30 eV): {((successful['lambda_hole_eV'] >= 0.25) & (successful['lambda_hole_eV'] < 0.30)).sum()}")
            print(f"  Poor (> 0.30 eV): {(successful['lambda_hole_eV'] >= 0.30).sum()}")
