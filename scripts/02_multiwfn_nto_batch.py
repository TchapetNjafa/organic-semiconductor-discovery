#!/usr/bin/env python3
"""
Multiwfn NTO Analysis - Tier 2
Automated NTO analysis for all successful molecules from Phase 2
"""

import os
import subprocess
import pandas as pd
from pathlib import Path
import time
import shutil
from joblib import Parallel, delayed
import multiprocessing

def prepare_wavefunction_xtb(xyz_file, work_dir):
    """Generate wavefunction file using xTB"""
    mol_name = Path(xyz_file).stem
    mol_dir = Path(work_dir) / mol_name
    mol_dir.mkdir(parents=True, exist_ok=True)
    
    # Copy XYZ to work directory
    work_xyz = mol_dir / 'mol.xyz'
    shutil.copy(xyz_file, work_xyz)
    
    try:
        # Run xtb4stda to generate wfn.xtb
        result = subprocess.run(
            ['xtb4stda', 'mol.xyz'],
            cwd=str(mol_dir),
            capture_output=True,
            text=True,
            timeout=300
        )
        
        wfn_file = mol_dir / 'wfn.xtb'
        if wfn_file.exists():
            return wfn_file
        return None
        
    except Exception as e:
        return None

def run_multiwfn_nto(wfn_file, mol_name, excited_state=1):
    """Run sTDA + Multiwfn NTO analysis"""
    mol_dir = wfn_file.parent
    
    try:
        # Step 1: Run sTDA to generate NTO molden file
        stda_result = subprocess.run(
            ['/home/tchapet/stda_v1-6/stda', '-xtb', 'wfn.xtb', '-nto', str(excited_state)],
            cwd=str(mol_dir),
            capture_output=True,
            text=True,
            timeout=300
        )
        
        # Check for NTO molden file
        nto_molden = mol_dir / f'nto{excited_state:03d}.molden'
        if not nto_molden.exists():
            return {'mol_id': mol_name, 'status': 'nto_generation_failed'}
        
        # Step 2: Run Multiwfn NTO analysis
        multiwfn_input = f"""18
1
{excited_state}
1
0
q
"""
        
        multiwfn_result = subprocess.run(
            ['Multiwfn', str(nto_molden)],
            input=multiwfn_input,
            capture_output=True,
            text=True,
            cwd=str(mol_dir),
            timeout=300
        )
        
        # Save Multiwfn output
        output_file = mol_dir / 'multiwfn_nto.log'
        with open(output_file, 'w') as f:
            f.write(multiwfn_result.stdout)
        
        # Parse results
        return parse_multiwfn_nto_output(multiwfn_result.stdout, mol_name)
        
    except subprocess.TimeoutExpired:
        return {'mol_id': mol_name, 'status': 'timeout'}
    except Exception as e:
        return {'mol_id': mol_name, 'status': f'error: {str(e)}'}

def parse_multiwfn_nto_output(output, mol_name):
    """Extract NTO properties from Multiwfn output"""
    data = {'mol_id': mol_name, 'status': 'success'}
    
    lines = output.split('\n')
    for i, line in enumerate(lines):
        # Lambda values (NTO weights)
        if 'Lambda' in line and 'Pair' in line:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    pair_num = int(parts[1])
                    lambda_val = float(parts[3])
                    data[f'Lambda_{pair_num}'] = lambda_val
                except:
                    pass
        
        # Sr index (electron-hole overlap)
        if 'Sr index' in line or 'Overlap integral' in line:
            parts = line.split()
            try:
                data['Sr_index'] = float(parts[-1])
            except:
                pass
        
        # D index (charge-transfer distance)
        if 'D index' in line or 'CT distance' in line:
            parts = line.split()
            try:
                # Find the value in Angstrom
                for j, p in enumerate(parts):
                    if 'Angstrom' in p or 'A' == p:
                        data['D_index_Angstrom'] = float(parts[j-1])
                        break
            except:
                pass
        
        # Transition dipole moment
        if 'Transition dipole moment' in line:
            try:
                if i+1 < len(lines):
                    parts = lines[i+1].split()
                    if len(parts) >= 4:
                        data['mu_total_Debye'] = float(parts[-1])
            except:
                pass
    
    return data

def process_single_molecule(xyz_file, work_dir, cleanup=True):
    """Process single molecule: xtb4stda + sTDA + Multiwfn"""
    mol_name = Path(xyz_file).stem
    
    try:
        # Generate wavefunction
        wfn_file = prepare_wavefunction_xtb(xyz_file, work_dir)
        if not wfn_file:
            return {'mol_id': mol_name, 'status': 'wfn_failed'}
        
        # Run NTO analysis
        nto_data = run_multiwfn_nto(wfn_file, mol_name, excited_state=1)
        
        # Cleanup if requested
        if cleanup and nto_data.get('status') == 'success':
            cleanup_temp_files(wfn_file.parent)
        
        return nto_data
        
    except Exception as e:
        return {'mol_id': mol_name, 'status': f'error: {str(e)}'}

def cleanup_temp_files(mol_dir):
    """Remove temporary files, keep only essential results"""
    mol_dir = Path(mol_dir)
    
    # Files to remove
    temp_patterns = [
        'wfn.xtb',
        'charges',
        'wbo',
        'xtbrestart',
        'xtbtopo.mol',
        'tda.dat',
        'nto*.molden'
    ]
    
    for pattern in temp_patterns:
        for file in mol_dir.glob(pattern):
            try:
                file.unlink()
            except:
                pass

def batch_nto_analysis(xyz_dir, output_csv, n_jobs=-1, cleanup=True):
    """Process all molecules for NTO analysis in parallel"""
    xyz_files = list(Path(xyz_dir).glob('*.xyz'))
    
    print(f"Found {len(xyz_files)} molecules to process")
    print(f"Using {n_jobs if n_jobs > 0 else multiprocessing.cpu_count()} CPU cores")
    print(f"Cleanup: {'Enabled' if cleanup else 'Disabled'}")
    
    work_dir = Path(__file__).parent.parent / 'data' / 'multiwfn_results'
    work_dir.mkdir(parents=True, exist_ok=True)
    
    start_time = time.time()
    
    # Parallel processing
    results = Parallel(n_jobs=n_jobs, verbose=10)(
        delayed(process_single_molecule)(str(xyz_file), work_dir, cleanup)
        for xyz_file in xyz_files
    )
    
    # Save final results
    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    print(f"\nCompleted! Results saved to {output_csv}")
    print(f"Total time: {(time.time() - start_time)/3600:.2f} hours")
    
    return df

if __name__ == '__main__':
    # Configuration
    XYZ_DIR = '../data/input_geometries'  # All molecules
    OUTPUT_CSV = '../analysis/multiwfn_nto_analysis.csv'
    
    print("="*60)
    print("Phase 3: Multiwfn NTO Analysis")
    print("="*60)
    print("")
    
    # Run batch NTO analysis
    # n_jobs: -1 for all CPUs, or specify number (e.g., 4)
    # cleanup: True to save disk space, False to keep all files
    df = batch_nto_analysis(XYZ_DIR, OUTPUT_CSV, n_jobs=-1, cleanup=True)
    
    # Statistics
    if len(df) > 0:
        print("\n=== NTO Analysis Summary ===")
        if 'status' in df.columns:
            print(f"Successful: {(df['status'] == 'success').sum()}")
            print(f"Failed: {(df['status'] != 'success').sum()}")
        if 'Sr_index' in df.columns:
            successful = df[df['status'] == 'success']
            if len(successful) > 0:
                print(f"\nSr index: {successful['Sr_index'].mean():.3f} ± {successful['Sr_index'].std():.3f}")
                print(f"  Range: {successful['Sr_index'].min():.3f} - {successful['Sr_index'].max():.3f}")
        if 'D_index_Angstrom' in df.columns:
            successful = df[df['status'] == 'success']
            if len(successful) > 0:
                print(f"\nD index: {successful['D_index_Angstrom'].mean():.2f} ± {successful['D_index_Angstrom'].std():.2f} Å")
                print(f"  Range: {successful['D_index_Angstrom'].min():.2f} - {successful['D_index_Angstrom'].max():.2f} Å")
        if 'Lambda_1' in df.columns:
            successful = df[df['status'] == 'success']
            if len(successful) > 0:
                print(f"\nLambda_1: {successful['Lambda_1'].mean():.3f} ± {successful['Lambda_1'].std():.3f}")
                print(f"  Single-excitation character (>0.8): {(successful['Lambda_1'] > 0.8).sum()}")
    else:
        print("\nNo molecules processed.")
