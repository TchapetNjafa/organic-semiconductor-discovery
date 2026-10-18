#!/usr/bin/env python3
"""
PySCF TD-DFT Validation - Tier 3
High-accuracy excited-state calculations with checkpoint support
Note: PySCF is memory-intensive, use limited parallelization
"""

import os
import numpy as np
import pandas as pd
from pyscf import gto, dft, tddft, lib
from pathlib import Path
from joblib import Parallel, delayed
import multiprocessing

# Set temp directory to current work directory
WORK_DIR = Path(__file__).parent.parent / 'data' / 'pyscf_results'
WORK_DIR.mkdir(parents=True, exist_ok=True)
os.environ['PYSCF_TMPDIR'] = str(WORK_DIR / 'tmp')
lib.param.TMPDIR = os.environ['PYSCF_TMPDIR']
Path(lib.param.TMPDIR).mkdir(exist_ok=True)

def run_tddft_pyscf(xyz_file, functional='cam-b3lyp', basis='def2-svp', nstates=5, cleanup=True):
    """
    Run TD-DFT calculation with PySCF and checkpoint support
    
    Parameters:
    -----------
    cleanup : bool
        Remove temporary files after calculation (default=True)
    """
    mol_name = Path(xyz_file).stem
    mol_dir = WORK_DIR / mol_name
    mol_dir.mkdir(parents=True, exist_ok=True)
    
    chkfile = mol_dir / f"{mol_name}.chk"
    
    # Read geometry
    mol = gto.M(atom=str(xyz_file), basis=basis, verbose=0)
    
    # Ground-state DFT with checkpoint
    mf = dft.RKS(mol)
    mf.xc = functional
    mf.chkfile = str(chkfile)
    
    # Try to load checkpoint
    if chkfile.exists():
        try:
            mf.__dict__.update(lib.chkfile.load(str(chkfile), 'scf'))
            print(f"  Loaded checkpoint for {mol_name}")
        except:
            mf.kernel()
    else:
        mf.kernel()
    
    if not mf.converged:
        return {'mol_id': mol_name, 'status': 'scf_failed'}
    
    # TD-DFT
    td = tddft.TDDFT(mf)
    td.nstates = nstates
    
    try:
        td.kernel()
    except Exception as e:
        # Try with TDA (Tamm-Dancoff approximation) if full TD-DFT fails
        print(f"  TD-DFT failed for {mol_name}, trying TDA...")
        try:
            td = tddft.TDA(mf)
            td.nstates = nstates
            td.kernel()
        except Exception as e2:
            return {'mol_id': mol_name, 'status': f'tddft_failed: {str(e2)[:50]}'}
    
    # Extract results
    data = {
        'mol_id': mol_name,
        'status': 'success',
        'functional': functional,
        'basis': basis,
        'E_HOMO': mf.mo_energy[mf.mo_occ > 0][-1] * 27.2114,  # to eV
        'E_LUMO': mf.mo_energy[mf.mo_occ == 0][0] * 27.2114,
        'E_gap_KS': (mf.mo_energy[mf.mo_occ == 0][0] - mf.mo_energy[mf.mo_occ > 0][-1]) * 27.2114
    }
    
    # Excitation energies and oscillator strengths
    for i in range(min(nstates, len(td.e))):
        data[f'S{i+1}_eV'] = td.e[i] * 27.2114
        data[f'S{i+1}_f'] = td.oscillator_strength()[i]
    
    # Analyze NTO for S1
    try:
        if len(td.xy) > 0:
            nto_data = analyze_nto_pyscf(td, state=0)
            data.update(nto_data)
    except Exception as e:
        print(f"  NTO analysis failed for {mol_name}: {str(e)[:50]}")
    
    # Cleanup temporary files if requested
    if cleanup:
        cleanup_pyscf_temp(mol_dir)
    
    return data

def cleanup_pyscf_temp(mol_dir):
    """
    Remove PySCF temporary files, keep checkpoint
    
    Keeps:
    - *.chk (checkpoint for restart)
    
    Removes:
    - tmp/ directory (scratch files)
    """
    import shutil
    mol_dir = Path(mol_dir)
    tmp_dir = mol_dir / 'tmp'
    
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir, ignore_errors=True)

def analyze_nto_pyscf(td, state=0):
    """Analyze NTO for given excited state"""
    # Get transition density matrix
    X, Y = td.xy[state]
    
    # Compute NTO via SVD
    U, s, Vt = np.linalg.svd(X + Y)
    
    nto_data = {}
    # Lambda values (singular values squared)
    for i in range(min(5, len(s))):
        nto_data[f'Lambda_{i+1}'] = s[i]**2
    
    # Single-excitation character
    nto_data['single_exc_character'] = (s[0]**2).sum()
    
    return nto_data

def calculate_reorganization_energy(xyz_file, functional='b3lyp', basis='6-31g*'):
    """Calculate reorganization energies for charge transport"""
    mol_name = Path(xyz_file).stem
    
    # Neutral geometry
    mol_n = gto.M(atom=str(xyz_file), basis=basis, verbose=0)
    mf_n = dft.RKS(mol_n)
    mf_n.xc = functional
    E_n = mf_n.kernel()
    
    # Cation at neutral geometry
    mol_c = gto.M(atom=str(xyz_file), basis=basis, charge=1, spin=1, verbose=0)
    mf_c = dft.RKS(mol_c)
    mf_c.xc = functional
    E_c_n = mf_c.kernel()
    
    # Optimize cation geometry (simplified - single point here)
    E_c_c = E_c_n  # Would need geometry optimization
    
    # Reorganization energy
    lambda_hole = (E_c_n - E_n) - (E_c_c - E_n)
    
    return {
        'mol_id': mol_name,
        'lambda_hole_eV': lambda_hole * 27.2114
    }

def batch_pyscf_analysis(xyz_dir, output_csv, max_molecules=20, n_jobs=1, cleanup=True):
    """
    Process top molecules with PySCF
    
    Parameters:
    -----------
    n_jobs : int
        Number of parallel jobs (default=1, PySCF is memory-intensive)
        Recommended: 1-4 depending on available RAM
    cleanup : bool
        Remove temporary files after calculation (default=True)
    """
    xyz_files = list(Path(xyz_dir).glob('*.xyz'))[:max_molecules]
    
    print(f"Processing {len(xyz_files)} molecules with PySCF TD-DFT")
    print(f"Using {n_jobs} parallel job(s)")
    print(f"Cleanup: {'Enabled' if cleanup else 'Disabled'}")
    print("Note: PySCF is memory-intensive, use n_jobs=1-4 max")
    
    # Parallel processing with limited jobs
    results = Parallel(n_jobs=n_jobs, verbose=10)(
        delayed(run_tddft_pyscf)(xyz_file, functional='cam-b3lyp', basis='def2-svp', cleanup=cleanup)
        for xyz_file in xyz_files
    )
    
    # Save final results
    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    print(f"\nCompleted! Results saved to {output_csv}")
    
    return df

if __name__ == '__main__':
    # Configuration
    XYZ_DIR = '../data/input_geometries'  # Top 7-20 molecules
    OUTPUT_CSV = '../analysis/pyscf_detailed_analysis.csv'
    
    # Run PySCF analysis
    # n_jobs=1 recommended (memory-intensive), increase if you have >32GB RAM
    # Set cleanup=False to keep all temporary files for debugging
    df = batch_pyscf_analysis(XYZ_DIR, OUTPUT_CSV, max_molecules=7, n_jobs=1, cleanup=True)
    
    # Statistics
    if len(df) > 0:
        print("\n=== PySCF TD-DFT Summary ===")
        if 'status' in df.columns:
            print(f"Successful: {(df['status'] == 'success').sum()}")
        if 'S1_eV' in df.columns:
            print(f"S1 Energy: {df['S1_eV'].mean():.2f} ± {df['S1_eV'].std():.2f} eV")
            print(f"S1 f: {df['S1_f'].mean():.3f} ± {df['S1_f'].std():.3f}")
