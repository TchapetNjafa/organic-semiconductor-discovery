#!/usr/bin/env python3
"""
Phase 8: Intrinsic Bond Orbital (IBO) Analysis
Provides localized orbital analysis complementary to NBO
IBOs are chemically intuitive and computationally efficient
"""

import os
import pandas as pd
from pyscf import gto, scf, lo, lib
from pathlib import Path
import numpy as np

# Set temp directory
WORK_DIR = Path(__file__).parent.parent / 'data' / 'ibo_analysis'
WORK_DIR.mkdir(parents=True, exist_ok=True)
os.environ['PYSCF_TMPDIR'] = str(WORK_DIR / 'tmp')
lib.param.TMPDIR = os.environ['PYSCF_TMPDIR']
Path(lib.param.TMPDIR).mkdir(exist_ok=True)

def run_ibo_analysis(xyz_file, functional='b3lyp', basis='6-31g*'):
    """
    Perform IBO analysis using PySCF
    
    Returns:
    - Number of IBOs
    - IBO localization metrics
    - Orbital centers
    """
    mol_name = Path(xyz_file).stem
    
    try:
        # Read geometry
        with open(xyz_file, 'r') as f:
            lines = f.readlines()
        
        natoms = int(lines[0].strip())
        atom_str = ''.join(lines[2:2+natoms])
        
        # Build molecule
        mol = gto.M(atom=atom_str, basis=basis, verbose=0)
        
        # SCF calculation
        mf = scf.RHF(mol)
        mf.kernel()
        
        if not mf.converged:
            return {
                'mol_id': mol_name,
                'status': 'scf_failed',
                'n_atoms': natoms,
                'n_basis': mol.nao
            }
        
        # IBO localization (occupied orbitals only)
        n_occ = int(mol.nelectron // 2)
        ibo_orbs = lo.ibo.ibo(mol, mf.mo_coeff[:, :n_occ])
        
        # Calculate IBO centers (expectation value of position)
        ibo_centers = []
        for i in range(ibo_orbs.shape[1]):
            orb = ibo_orbs[:, i]
            # Mulliken population for this IBO
            pop = mf.mulliken_pop(mol, np.outer(orb, orb), s=mf.get_ovlp())[1]
            # Weighted center
            coords = mol.atom_coords()
            center = np.average(coords, axis=0, weights=pop)
            ibo_centers.append(center)
        
        ibo_centers = np.array(ibo_centers)
        
        # Calculate localization metric (spread of IBO centers)
        center_spread = np.std(np.linalg.norm(ibo_centers, axis=1))
        
        # Calculate average IBO-atom distance (measure of localization)
        atom_coords = mol.atom_coords()
        min_distances = []
        for center in ibo_centers:
            dists = np.linalg.norm(atom_coords - center, axis=1)
            min_distances.append(np.min(dists))
        avg_ibo_atom_dist = np.mean(min_distances)
        
        return {
            'mol_id': mol_name,
            'status': 'success',
            'n_atoms': natoms,
            'n_basis': mol.nao,
            'n_electrons': mol.nelectron,
            'n_ibos': ibo_orbs.shape[1],
            'E_total_hartree': mf.e_tot,
            'ibo_center_spread_angstrom': center_spread,
            'avg_ibo_atom_distance_angstrom': avg_ibo_atom_dist,
            'localization_quality': 'Good' if avg_ibo_atom_dist < 0.5 else 'Moderate'
        }
        
    except Exception as e:
        return {
            'mol_id': mol_name,
            'status': f'error: {str(e)[:50]}'
        }

def batch_ibo_analysis(xyz_dir, output_csv, mol_list=None):
    """
    Run IBO analysis for specified molecules
    """
    if mol_list:
        xyz_files = [Path(xyz_dir) / f"{mol_id}.xyz" for mol_id in mol_list]
        xyz_files = [f for f in xyz_files if f.exists()]
    else:
        xyz_files = list(Path(xyz_dir).glob('*.xyz'))[:5]
    
    print(f"Processing {len(xyz_files)} molecules for IBO analysis")
    print(f"Estimated time: ~{len(xyz_files) * 1.3:.0f} minutes\n")
    
    results = []
    for i, xyz_file in enumerate(xyz_files, 1):
        print(f"[{i}/{len(xyz_files)}] Analyzing {xyz_file.stem}...")
        result = run_ibo_analysis(str(xyz_file))
        results.append(result)
        
        if result['status'] == 'success':
            print(f"  ✓ {result['n_ibos']} IBOs, localization: {result['localization_quality']}")
        else:
            print(f"  ✗ {result['status']}")
    
    # Save results
    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    print(f"\n{'='*60}")
    print(f"Completed! Results saved to {output_csv}")
    print(f"{'='*60}")
    
    return df

if __name__ == '__main__':
    XYZ_DIR = '../data/input_geometries'
    OUTPUT_CSV = '../analysis/ibo_analysis.csv'
    
    # Top molecules from manuscript
    KEY_MOLECULES = [17851, 19531, 14380, 20778, 977, 7801, 18985]
    
    print("="*60)
    print("Phase 8: Intrinsic Bond Orbital (IBO) Analysis")
    print("="*60)
    print("")
    print("IBO provides:")
    print("  - Localized orbital representation")
    print("  - Chemically intuitive bonding picture")
    print("  - Complementary to NBO analysis")
    print("")
    
    # Run analysis
    df = batch_ibo_analysis(XYZ_DIR, OUTPUT_CSV, mol_list=KEY_MOLECULES)
    
    # Statistics
    if len(df) > 0:
        print("\n=== IBO Analysis Summary ===")
        successful = df[df['status'] == 'success']
        if len(successful) > 0:
            print(f"Successful: {len(successful)}/{len(df)}")
            print(f"Average IBOs per molecule: {successful['n_ibos'].mean():.1f}")
            print(f"Average localization distance: {successful['avg_ibo_atom_distance_angstrom'].mean():.3f} Å")
            print(f"\nTop 3 most localized:")
            top3 = successful.nsmallest(3, 'avg_ibo_atom_distance_angstrom')
            for _, row in top3.iterrows():
                print(f"  Mol {row['mol_id']}: {row['avg_ibo_atom_distance_angstrom']:.3f} Å")
