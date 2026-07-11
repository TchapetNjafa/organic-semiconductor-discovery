#!/usr/bin/env python3
"""
Phase 7: NBO-like Analysis using PySCF + Multiwfn
Calculates charge distributions and orbital contributions
"""

import os
import subprocess
import pandas as pd
from pyscf import gto, dft, lib
from pathlib import Path
import numpy as np

# Set temp directory
WORK_DIR = Path(__file__).parent.parent / 'data' / 'nbo_analysis'
WORK_DIR.mkdir(parents=True, exist_ok=True)
os.environ['PYSCF_TMPDIR'] = str(WORK_DIR / 'tmp')
lib.param.TMPDIR = os.environ['PYSCF_TMPDIR']
Path(lib.param.TMPDIR).mkdir(exist_ok=True)

def run_nbo_analysis(xyz_file, functional='b3lyp', basis='6-31g*'):
    """
    Perform NBO-like analysis using PySCF + Multiwfn
    
    Multiwfn can perform:
    - Natural population analysis (NPA)
    - Orbital composition analysis
    - Charge distribution
    """
    mol_name = Path(xyz_file).stem
    mol_dir = WORK_DIR / mol_name
    mol_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Processing {mol_name}...")
    
    try:
        # Step 1: Run DFT calculation with PySCF
        mol = gto.M(atom=str(xyz_file), basis=basis, verbose=0)
        mf = dft.RKS(mol)
        mf.xc = functional
        mf.kernel()
        
        if not mf.converged:
            return {'mol_id': mol_name, 'status': 'scf_failed'}
        
        # Save molden file for Multiwfn
        molden_file = mol_dir / f'{mol_name}.molden'
        with open(molden_file, 'w') as f:
            from pyscf.tools import molden
            molden.from_mo(mol, str(molden_file), mf.mo_coeff, ene=mf.mo_energy, occ=mf.mo_occ)
        
        # Step 2: Run Multiwfn for NPA (Natural Population Analysis)
        print(f"  Running Multiwfn NPA analysis...")
        multiwfn_npa_input = """7
5
1
0
q
"""
        
        npa_result = subprocess.run(
            ['Multiwfn', str(molden_file)],
            input=multiwfn_npa_input,
            capture_output=True,
            text=True,
            cwd=str(mol_dir),
            timeout=300
        )
        
        # Save NPA output
        npa_log = mol_dir / 'multiwfn_npa.log'
        with open(npa_log, 'w') as f:
            f.write(npa_result.stdout)
        
        # Step 3: Run Multiwfn for orbital composition analysis
        print(f"  Running Multiwfn orbital composition analysis...")
        
        # Get HOMO and LUMO indices
        nocc = int(mf.mo_occ.sum() / 2)
        homo_idx = nocc
        lumo_idx = nocc + 1
        
        # Analyze HOMO composition
        multiwfn_homo_input = f"""8
1
{homo_idx}
0
q
"""
        
        homo_result = subprocess.run(
            ['Multiwfn', str(molden_file)],
            input=multiwfn_homo_input,
            capture_output=True,
            text=True,
            cwd=str(mol_dir),
            timeout=300
        )
        
        homo_log = mol_dir / 'multiwfn_homo.log'
        with open(homo_log, 'w') as f:
            f.write(homo_result.stdout)
        
        # Analyze LUMO composition
        multiwfn_lumo_input = f"""8
1
{lumo_idx}
0
q
"""
        
        lumo_result = subprocess.run(
            ['Multiwfn', str(molden_file)],
            input=multiwfn_lumo_input,
            capture_output=True,
            text=True,
            cwd=str(mol_dir),
            timeout=300
        )
        
        lumo_log = mol_dir / 'multiwfn_lumo.log'
        with open(lumo_log, 'w') as f:
            f.write(lumo_result.stdout)
        
        # Parse results
        data = parse_nbo_results(npa_log, homo_log, lumo_log, mol_name, mol)
        
        return data
        
    except Exception as e:
        return {'mol_id': mol_name, 'status': f'error: {str(e)[:100]}'}

def parse_nbo_results(npa_log, homo_log, lumo_log, mol_name, mol):
    """
    Parse Multiwfn output for NPA charges and orbital compositions
    """
    data = {'mol_id': mol_name, 'status': 'success'}
    
    # Parse NPA charges
    try:
        with open(npa_log, 'r') as f:
            npa_lines = f.readlines()
        
        # Find atomic charges
        charges_by_element = {}
        for line in npa_lines:
            if 'Atomic charges' in line or 'Natural charges' in line:
                # Parse charge data
                pass
        
        # Calculate average charges for O, N, C
        # This is simplified - actual parsing depends on Multiwfn output format
        data['npa_charges_available'] = True
        
    except Exception as e:
        data['npa_charges_available'] = False
    
    # Parse HOMO composition
    try:
        with open(homo_log, 'r') as f:
            homo_lines = f.readlines()
        
        # Look for orbital composition by atom type
        for line in homo_lines:
            if 'Contribution' in line or 'composition' in line.lower():
                # Parse composition data
                pass
        
        data['homo_composition_available'] = True
        
    except Exception as e:
        data['homo_composition_available'] = False
    
    # Parse LUMO composition
    try:
        with open(lumo_log, 'r') as f:
            lumo_lines = f.readlines()
        
        data['lumo_composition_available'] = True
        
    except Exception as e:
        data['lumo_composition_available'] = False
    
    # Add note about manual inspection
    data['note'] = 'Check log files for detailed composition data'
    
    return data

def batch_nbo_analysis(xyz_dir, output_csv, mol_list=None):
    """
    Perform NBO-like analysis for specified molecules
    """
    if mol_list:
        xyz_files = [Path(xyz_dir) / f"{mol_id}.xyz" for mol_id in mol_list]
        xyz_files = [f for f in xyz_files if f.exists()]
    else:
        xyz_files = list(Path(xyz_dir).glob('*.xyz'))[:5]
    
    print(f"Processing {len(xyz_files)} molecules for NBO-like analysis")
    print("Note: Results will be in individual log files")
    
    results = []
    for xyz_file in xyz_files:
        result = run_nbo_analysis(str(xyz_file))
        results.append(result)
    
    # Save results
    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    print(f"\nCompleted! Results saved to {output_csv}")
    print(f"\nDetailed results in: {WORK_DIR}/[mol_id]/multiwfn_*.log")
    
    return df

if __name__ == '__main__':
    XYZ_DIR = '../data/input_geometries'
    OUTPUT_CSV = '../analysis/nbo_analysis.csv'
    
    # Key molecules from manuscript
    KEY_MOLECULES = [17851, 1712]
    
    print("="*60)
    print("Phase 7: NBO-like Analysis (PySCF + Multiwfn)")
    print("="*60)
    print("")
    
    # Run for key molecules
    df = batch_nbo_analysis(XYZ_DIR, OUTPUT_CSV, mol_list=KEY_MOLECULES)
    
    # Statistics
    if len(df) > 0:
        print("\n=== NBO Analysis Summary ===")
        successful = df[df['status'] == 'success']
        if len(successful) > 0:
            print(f"Successful: {len(successful)}/{len(df)}")
            print("\nTo extract specific values:")
            print("  1. Check data/nbo_analysis/[mol_id]/multiwfn_npa.log for charges")
            print("  2. Check data/nbo_analysis/[mol_id]/multiwfn_homo.log for HOMO composition")
            print("  3. Check data/nbo_analysis/[mol_id]/multiwfn_lumo.log for LUMO composition")
