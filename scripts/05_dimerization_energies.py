#!/usr/bin/env python3
"""
Phase 5: Dimerization Energy Calculations
Calculate aggregation tendency for solid-state morphology assessment
"""

import os
import pandas as pd
from pyscf import gto, dft, lib
from pathlib import Path
import numpy as np

# Set temp directory
WORK_DIR = Path(__file__).parent.parent / 'data' / 'dimerization_results'
WORK_DIR.mkdir(parents=True, exist_ok=True)
os.environ['PYSCF_TMPDIR'] = str(WORK_DIR / 'tmp')
lib.param.TMPDIR = os.environ['PYSCF_TMPDIR']
Path(lib.param.TMPDIR).mkdir(exist_ok=True)

def calculate_dimerization_energy(xyz_file, functional='b3lyp', basis='6-31g*', separation=3.5):
    """
    Calculate dimerization energy for π-π stacking
    
    E_dimer = E(dimer) - 2*E(monomer)
    
    Parameters:
    -----------
    separation : float
        Initial π-π stacking distance in Angstrom (typical: 3.3-3.5)
    """
    mol_name = Path(xyz_file).stem
    
    print(f"Processing {mol_name}...")
    
    try:
        # Read monomer geometry
        with open(xyz_file, 'r') as f:
            lines = f.readlines()
        
        natoms = int(lines[0].strip())
        atoms = []
        coords = []
        
        for line in lines[2:2+natoms]:
            parts = line.split()
            if len(parts) >= 4:
                atoms.append(parts[0])
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
        
        coords = np.array(coords)
        
        # Calculate monomer energy
        atom_str = '\n'.join([f"{atom} {x:.6f} {y:.6f} {z:.6f}" 
                              for atom, (x, y, z) in zip(atoms, coords)])
        
        mol_monomer = gto.M(atom=atom_str, basis=basis, verbose=0)
        mf_monomer = dft.RKS(mol_monomer)
        mf_monomer.xc = functional
        E_monomer = mf_monomer.kernel()
        
        if not mf_monomer.converged:
            return {'mol_id': mol_name, 'status': 'monomer_scf_failed'}
        
        # Create dimer by stacking along z-axis
        # Shift second molecule by separation distance
        coords_dimer = np.vstack([coords, coords + np.array([0, 0, separation])])
        atoms_dimer = atoms * 2
        
        atom_str_dimer = '\n'.join([f"{atom} {x:.6f} {y:.6f} {z:.6f}" 
                                     for atom, (x, y, z) in zip(atoms_dimer, coords_dimer)])
        
        mol_dimer = gto.M(atom=atom_str_dimer, basis=basis, verbose=0)
        mf_dimer = dft.RKS(mol_dimer)
        mf_dimer.xc = functional
        E_dimer = mf_dimer.kernel()
        
        if not mf_dimer.converged:
            return {'mol_id': mol_name, 'status': 'dimer_scf_failed'}
        
        # Calculate dimerization energy (in kcal/mol)
        E_dim = (E_dimer - 2 * E_monomer) * 627.509  # Hartree to kcal/mol
        
        # Classify aggregation tendency
        if E_dim < -6.0:
            classification = "Strong aggregation"
        elif E_dim < -3.0:
            classification = "Moderate aggregation"
        elif E_dim < -1.0:
            classification = "Mild aggregation"
        else:
            classification = "Weak/no aggregation"
        
        return {
            'mol_id': mol_name,
            'status': 'success',
            'E_dimerization_kcal_mol': E_dim,
            'separation_angstrom': separation,
            'aggregation_class': classification,
            'E_monomer_hartree': E_monomer,
            'E_dimer_hartree': E_dimer
        }
        
    except Exception as e:
        return {'mol_id': mol_name, 'status': f'error: {str(e)[:100]}'}

def batch_dimerization_analysis(xyz_dir, output_csv, mol_list=None):
    """
    Calculate dimerization energies for specified molecules
    
    Parameters:
    -----------
    mol_list : list
        List of molecule IDs (e.g., [17851, 7801])
    """
    if mol_list:
        xyz_files = [Path(xyz_dir) / f"{mol_id}.xyz" for mol_id in mol_list]
        xyz_files = [f for f in xyz_files if f.exists()]
    else:
        xyz_files = list(Path(xyz_dir).glob('*.xyz'))[:10]  # Limit to 10 if no list
    
    print(f"Processing {len(xyz_files)} molecules for dimerization energy")
    print("Note: Sequential processing (memory intensive)")
    
    results = []
    for xyz_file in xyz_files:
        result = calculate_dimerization_energy(str(xyz_file))
        results.append(result)
    
    # Save results
    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    print(f"\nCompleted! Results saved to {output_csv}")
    
    return df

if __name__ == '__main__':
    XYZ_DIR = '../data/input_geometries'
    OUTPUT_CSV = '../analysis/dimerization_energies.csv'
    
    # Key molecules mentioned in manuscript
    KEY_MOLECULES = [17851, 7801]
    
    print("="*60)
    print("Phase 5: Dimerization Energy Calculations")
    print("="*60)
    print("")
    
    # Run for key molecules
    df = batch_dimerization_analysis(XYZ_DIR, OUTPUT_CSV, mol_list=KEY_MOLECULES)
    
    # Statistics
    if len(df) > 0:
        print("\n=== Dimerization Energy Summary ===")
        successful = df[df['status'] == 'success']
        if len(successful) > 0:
            print(f"Successful: {len(successful)}/{len(df)}")
            print(f"\nDimerization energies:")
            for _, row in successful.iterrows():
                print(f"  Mol {row['mol_id']}: {row['E_dimerization_kcal_mol']:.2f} kcal/mol ({row['aggregation_class']})")
            
            print(f"\nMean: {successful['E_dimerization_kcal_mol'].mean():.2f} kcal/mol")
            print(f"Range: {successful['E_dimerization_kcal_mol'].min():.2f} to {successful['E_dimerization_kcal_mol'].max():.2f} kcal/mol")
