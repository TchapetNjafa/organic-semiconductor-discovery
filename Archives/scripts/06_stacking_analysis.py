#!/usr/bin/env python3
"""
Phase 6: Simple π-π Stacking Analysis
Alternative to full crystal structure prediction
Analyzes optimal stacking distance and slip angle for molecular dimers
"""

import os
import pandas as pd
from pyscf import gto, dft, lib
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar

# Set temp directory
WORK_DIR = Path(__file__).parent.parent / 'data' / 'stacking_analysis'
WORK_DIR.mkdir(parents=True, exist_ok=True)
os.environ['PYSCF_TMPDIR'] = str(WORK_DIR / 'tmp')
lib.param.TMPDIR = os.environ['PYSCF_TMPDIR']
Path(lib.param.TMPDIR).mkdir(exist_ok=True)

def calculate_stacking_energy(xyz_file, separation, slip_angle=0, functional='b3lyp', basis='6-31g*'):
    """
    Calculate π-π stacking energy at given separation and slip angle
    
    Parameters:
    -----------
    separation : float
        π-π stacking distance in Angstrom
    slip_angle : float
        Slip angle in degrees (0 = perfect stacking)
    """
    mol_name = Path(xyz_file).stem
    
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
        
        # Center molecule
        center = coords.mean(axis=0)
        coords_centered = coords - center
        
        # Calculate monomer energy
        atom_str = '\n'.join([f"{atom} {x:.6f} {y:.6f} {z:.6f}" 
                              for atom, (x, y, z) in zip(atoms, coords_centered)])
        
        mol_monomer = gto.M(atom=atom_str, basis=basis, verbose=0)
        mf_monomer = dft.RKS(mol_monomer)
        mf_monomer.xc = functional
        E_monomer = mf_monomer.kernel()
        
        if not mf_monomer.converged:
            return None, None
        
        # Create dimer with slip
        slip_x = separation * np.tan(np.radians(slip_angle))
        coords_dimer = np.vstack([
            coords_centered,
            coords_centered + np.array([slip_x, 0, separation])
        ])
        atoms_dimer = atoms * 2
        
        atom_str_dimer = '\n'.join([f"{atom} {x:.6f} {y:.6f} {z:.6f}" 
                                     for atom, (x, y, z) in zip(atoms_dimer, coords_dimer)])
        
        mol_dimer = gto.M(atom=atom_str_dimer, basis=basis, verbose=0)
        mf_dimer = dft.RKS(mol_dimer)
        mf_dimer.xc = functional
        E_dimer = mf_dimer.kernel()
        
        if not mf_dimer.converged:
            return None, None
        
        # Calculate interaction energy (in kcal/mol)
        E_int = (E_dimer - 2 * E_monomer) * 627.509
        
        return E_int, E_dimer
        
    except Exception as e:
        return None, None

def optimize_stacking(xyz_file, functional='b3lyp', basis='6-31g*'):
    """
    Find optimal π-π stacking distance and slip angle
    """
    mol_name = Path(xyz_file).stem
    
    print(f"Optimizing stacking for {mol_name}...")
    
    # Scan separation distances
    separations = np.linspace(3.0, 4.0, 11)  # 3.0 to 4.0 Å
    energies = []
    
    print("  Scanning separation distances...")
    for sep in separations:
        E_int, _ = calculate_stacking_energy(xyz_file, sep, slip_angle=0)
        if E_int is not None:
            energies.append(E_int)
            print(f"    {sep:.2f} Å: {E_int:.2f} kcal/mol")
        else:
            energies.append(np.nan)
    
    # Find optimal separation
    valid_idx = ~np.isnan(energies)
    if not any(valid_idx):
        return {'mol_id': mol_name, 'status': 'all_failed'}
    
    optimal_idx = np.nanargmin(energies)
    optimal_sep = separations[optimal_idx]
    optimal_energy = energies[optimal_idx]
    
    print(f"  Optimal separation: {optimal_sep:.2f} Å ({optimal_energy:.2f} kcal/mol)")
    
    # Scan slip angles at optimal separation
    print("  Scanning slip angles...")
    slip_angles = [0, 10, 20, 30]
    slip_energies = []
    
    for angle in slip_angles:
        E_int, _ = calculate_stacking_energy(xyz_file, optimal_sep, slip_angle=angle)
        if E_int is not None:
            slip_energies.append(E_int)
            print(f"    {angle}°: {E_int:.2f} kcal/mol")
        else:
            slip_energies.append(np.nan)
    
    # Find optimal slip angle
    optimal_slip_idx = np.nanargmin(slip_energies)
    optimal_slip = slip_angles[optimal_slip_idx]
    optimal_slip_energy = slip_energies[optimal_slip_idx]
    
    print(f"  Optimal slip angle: {optimal_slip}° ({optimal_slip_energy:.2f} kcal/mol)")
    
    return {
        'mol_id': mol_name,
        'status': 'success',
        'optimal_separation_angstrom': optimal_sep,
        'optimal_slip_angle_degrees': optimal_slip,
        'interaction_energy_kcal_mol': optimal_slip_energy,
        'stacking_quality': 'Favorable' if optimal_sep >= 3.3 and optimal_sep <= 3.5 else 'Suboptimal'
    }

def batch_stacking_analysis(xyz_dir, output_csv, mol_list=None):
    """
    Analyze π-π stacking for specified molecules
    """
    if mol_list:
        xyz_files = [Path(xyz_dir) / f"{mol_id}.xyz" for mol_id in mol_list]
        xyz_files = [f for f in xyz_files if f.exists()]
    else:
        xyz_files = list(Path(xyz_dir).glob('*.xyz'))[:5]
    
    print(f"Processing {len(xyz_files)} molecules for stacking analysis")
    print("Note: This is computationally intensive")
    
    results = []
    for xyz_file in xyz_files:
        result = optimize_stacking(str(xyz_file))
        results.append(result)
    
    # Save results
    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    print(f"\nCompleted! Results saved to {output_csv}")
    
    return df

if __name__ == '__main__':
    XYZ_DIR = '../data/input_geometries'
    OUTPUT_CSV = '../analysis/stacking_analysis.csv'
    
    # Key molecule from manuscript
    KEY_MOLECULES = [17851]
    
    print("="*60)
    print("Phase 6: π-π Stacking Analysis")
    print("="*60)
    print("")
    
    # Run for key molecule
    df = batch_stacking_analysis(XYZ_DIR, OUTPUT_CSV, mol_list=KEY_MOLECULES)
    
    # Statistics
    if len(df) > 0:
        print("\n=== Stacking Analysis Summary ===")
        successful = df[df['status'] == 'success']
        if len(successful) > 0:
            print(f"Successful: {len(successful)}/{len(df)}")
            for _, row in successful.iterrows():
                print(f"\nMol {row['mol_id']}:")
                print(f"  Optimal separation: {row['optimal_separation_angstrom']:.2f} Å")
                print(f"  Optimal slip angle: {row['optimal_slip_angle_degrees']:.0f}°")
                print(f"  Interaction energy: {row['interaction_energy_kcal_mol']:.2f} kcal/mol")
                print(f"  Quality: {row['stacking_quality']}")
