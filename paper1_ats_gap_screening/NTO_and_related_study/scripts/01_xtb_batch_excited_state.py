#!/usr/bin/env python3
"""
xTB + sTDA Excited-State Batch Screening
Tier 1: High-throughput screening using xtb4stda + stda
With parallel processing support
"""

import os
import subprocess
import pandas as pd
from pathlib import Path
import time
import shutil
from joblib import Parallel, delayed
import multiprocessing

def run_xtb4stda_excited_state(xyz_file, output_dir, cleanup=True):
    """
    Run xtb4stda + sTDA for single molecule
    
    Parameters:
    -----------
    cleanup : bool
        If True, remove temporary files after extraction (default=True)
    """
    mol_name = Path(xyz_file).stem
    result_dir = Path(output_dir) / mol_name
    result_dir.mkdir(exist_ok=True)
    
    # Copy XYZ to work directory
    work_xyz = result_dir / 'mol.xyz'
    shutil.copy(xyz_file, work_xyz)
    
    try:
        # Run xtb4stda
        result = subprocess.run(
            ['xtb4stda', 'mol.xyz'],
            cwd=str(result_dir),
            capture_output=True,
            text=True,
            timeout=300
        )
        
        wfn_file = result_dir / 'wfn.xtb'
        if not wfn_file.exists():
            return {'mol_id': mol_name, 'status': 'wfn_failed'}
        
        # Run sTDA
        stda_result = subprocess.run(
            ['/home/tchapet/stda_v1-6/stda', '-xtb', 'wfn.xtb'],
            cwd=str(result_dir),
            capture_output=True,
            text=True,
            timeout=300
        )
        
        tda_file = result_dir / 'tda.dat'
        if not tda_file.exists():
            return {'mol_id': mol_name, 'status': 'stda_failed'}
        
        # Parse results
        data = parse_stda_output(tda_file, mol_name)
        
        # Cleanup temporary files if requested
        if cleanup and data.get('status') == 'success':
            cleanup_temp_files(result_dir)
        
        return data
        
    except subprocess.TimeoutExpired:
        return {'mol_id': mol_name, 'status': 'timeout'}
    except Exception as e:
        return {'mol_id': mol_name, 'status': f'error: {str(e)}'}

def parse_stda_output(tda_file, mol_name):
    """Extract excited-state properties from tda.dat"""
    data = {'mol_id': mol_name, 'status': 'success'}
    
    with open(tda_file, 'r') as f:
        lines = f.readlines()
    
    # Find excitation energies section
    for i, line in enumerate(lines):
        if 'excitation energies' in line.lower():
            # Parse next lines for S1-S10
            j = i + 3
            s_count = 1
            while j < len(lines) and s_count <= 10:
                parts = lines[j].split()
                if len(parts) >= 4 and parts[0].isdigit():
                    data[f'S{s_count}_eV'] = float(parts[1])
                    data[f'S{s_count}_nm'] = float(parts[2])
                    data[f'S{s_count}_f'] = float(parts[3])
                    s_count += 1
                j += 1
            break
    
    return data

def cleanup_temp_files(result_dir):
    """
    Remove temporary files, keep only essential results
    
    Keeps:
    - tda.dat (excitation energies)
    - mol.xyz (input geometry)
    
    Removes:
    - wfn.xtb (large wavefunction file)
    - xtb temporary files (charges, wbo, etc.)
    - log files
    """
    result_dir = Path(result_dir)
    
    # Files to remove
    temp_patterns = [
        'wfn.xtb',
        'charges',
        'wbo',
        'xtbrestart',
        'xtbtopo.mol',
        '*.log',
        'xtb.out',
        'xtb.trj'
    ]
    
    for pattern in temp_patterns:
        for file in result_dir.glob(pattern):
            try:
                file.unlink()
            except:
                pass

def batch_process(input_dir, output_dir, csv_output, n_jobs=-1, cleanup=True):
    """
    Process all XYZ files in input directory with parallel execution
    
    Parameters:
    -----------
    n_jobs : int
        Number of parallel jobs (-1 = all CPUs, 1 = sequential)
    cleanup : bool
        Remove temporary files after processing (default=True)
    """
    xyz_files = list(Path(input_dir).glob('*.xyz'))
    
    print(f"Found {len(xyz_files)} molecules to process")
    print(f"Using {n_jobs if n_jobs > 0 else multiprocessing.cpu_count()} CPU cores")
    print(f"Cleanup: {'Enabled' if cleanup else 'Disabled'}")
    
    start_time = time.time()
    
    # Parallel processing
    results = Parallel(n_jobs=n_jobs, verbose=10)(
        delayed(run_xtb4stda_excited_state)(str(xyz_file), output_dir, cleanup)
        for xyz_file in xyz_files
    )
    
    # Save final results
    df = pd.DataFrame(results)
    df.to_csv(csv_output, index=False)
    print(f"\nCompleted! Results saved to {csv_output}")
    print(f"Total time: {(time.time() - start_time)/3600:.2f} hours")
    
    return df

if __name__ == '__main__':
    # Configuration
    INPUT_DIR = '../data/input_geometries'
    OUTPUT_DIR = '../data/xtb_results'
    CSV_OUTPUT = '../analysis/xtb_excited_state_screening.csv'
    
    # Check if input directory exists and has files
    if not Path(INPUT_DIR).exists():
        print(f"Creating input directory: {INPUT_DIR}")
        Path(INPUT_DIR).mkdir(parents=True, exist_ok=True)
        print(f"\nPlease copy XYZ files to: {INPUT_DIR}")
        print("Example: cp /path/to/MVOTOOPV/figures/xTB_PRE_OPT_XYZ/*.xyz {INPUT_DIR}/")
        exit(1)
    
    xyz_files = list(Path(INPUT_DIR).glob('*.xyz'))
    if len(xyz_files) == 0:
        print(f"No XYZ files found in {INPUT_DIR}")
        print("\nTo populate with all molecules:")
        print(f"  cp '/home/tchapet/Post-Doc/ARTICLES DOSSIERS DES ARTICLES EN REDACTION/NOUVELS AXES DE RECHERCHE A REGARDER URGEMMENT/ARTICLES EN REDACTIONS/ARTICLE_MVOTO/MVOTOOPV/figures/xTB_PRE_OPT_XYZ'/*.xyz {INPUT_DIR}/")
        exit(1)
    
    # Run batch processing
    # Use n_jobs=1 for sequential, n_jobs=-1 for all CPUs, or specify number
    # Set cleanup=False to keep all temporary files for debugging
    df = batch_process(INPUT_DIR, OUTPUT_DIR, CSV_OUTPUT, n_jobs=-1, cleanup=True)
    
    # Quick statistics
    if len(df) > 0:
        print("\n=== Summary Statistics ===")
        if 'status' in df.columns:
            print(f"Successful: {(df['status'] == 'success').sum()}")
            print(f"Failed: {(df['status'] != 'success').sum()}")
        if 'S1_eV' in df.columns:
            print(f"\nS1 Energy: {df['S1_eV'].mean():.2f} ± {df['S1_eV'].std():.2f} eV")
            print(f"S1 f: {df['S1_f'].mean():.3f} ± {df['S1_f'].std():.3f}")
    else:
        print("\nNo molecules processed.")
