#!/usr/bin/env python3
"""
Test script to validate xTB, Multiwfn, and PySCF setup
Run this before starting the full workflow
"""

import subprocess
import sys
from pathlib import Path
import tempfile

def test_xtb():
    """Test xTB installation and excited-state capability"""
    print("Testing xTB...")
    
    # Create test molecule (benzene)
    test_xyz = """6
Benzene test
C        0.00000        1.40272        0.00000
C        1.21479        0.70136        0.00000
C        1.21479       -0.70136        0.00000
C        0.00000       -1.40272        0.00000
C       -1.21479       -0.70136        0.00000
C       -1.21479        0.70136        0.00000
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
        f.write(test_xyz)
        xyz_file = f.name
    
    try:
        # Test basic xTB
        result = subprocess.run(['xtb', '--version'], capture_output=True, text=True)
        print(f"  ✓ xTB version: {result.stdout.strip()}")
        
        # Test excited-state calculation
        result = subprocess.run(
            ['xtb', xyz_file, '--gfn', '2', '--stda'],
            capture_output=True, text=True, timeout=30
        )
        
        if 'excitation energies' in result.stdout.lower():
            print("  ✓ xTB excited-state calculation works")
            return True
        else:
            print("  ✗ xTB excited-state calculation failed")
            return False
            
    except FileNotFoundError:
        print("  ✗ xTB not found in PATH")
        return False
    except Exception as e:
        print(f"  ✗ xTB test failed: {e}")
        return False
    finally:
        Path(xyz_file).unlink(missing_ok=True)

def test_multiwfn():
    """Test Multiwfn installation"""
    print("\nTesting Multiwfn...")
    
    try:
        result = subprocess.run(['Multiwfn'], input='q\n', capture_output=True, text=True, timeout=5)
        if 'Multiwfn' in result.stdout or 'Tian Lu' in result.stdout:
            print("  ✓ Multiwfn is accessible")
            return True
        else:
            print("  ✗ Multiwfn response unexpected")
            return False
    except FileNotFoundError:
        print("  ✗ Multiwfn not found in PATH")
        return False
    except Exception as e:
        print(f"  ✗ Multiwfn test failed: {e}")
        return False

def test_pyscf():
    """Test PySCF installation and TD-DFT"""
    print("\nTesting PySCF...")
    
    try:
        from pyscf import gto, dft, tddft
        print("  ✓ PySCF imported successfully")
        
        # Quick test calculation
        mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g', verbose=0)
        mf = dft.RKS(mol)
        mf.xc = 'b3lyp'
        mf.kernel()
        
        td = tddft.TDDFT(mf)
        td.nstates = 3
        td.kernel()
        
        if td.converged:
            print("  ✓ PySCF TD-DFT calculation works")
            print(f"    S1 = {td.e[0]*27.2114:.2f} eV")
            return True
        else:
            print("  ✗ PySCF TD-DFT did not converge")
            return False
            
    except ImportError as e:
        print(f"  ✗ PySCF import failed: {e}")
        return False
    except Exception as e:
        print(f"  ✗ PySCF test failed: {e}")
        return False

def check_data_availability():
    """Check if PubChemQC data is available"""
    print("\nChecking data availability...")
    
    # Look for existing datasets
    data_paths = [
        '../MVOTOOPV/PUBCHEMQC_MULTIFONCTIONAL_2025/DATASET/dataset_pubchemqc_opv_17458.csv',
        '../MVOTOOPV/DATASET/',
    ]
    
    for path in data_paths:
        if Path(path).exists():
            print(f"  ✓ Found: {path}")
        else:
            print(f"  ✗ Not found: {path}")
    
    return True

def main():
    """Run all tests"""
    print("="*60)
    print("NTO Study Setup Validation")
    print("="*60)
    
    results = {
        'xTB': test_xtb(),
        'Multiwfn': test_multiwfn(),
        'PySCF': test_pyscf(),
        'Data': check_data_availability()
    }
    
    print("\n" + "="*60)
    print("Summary:")
    print("="*60)
    for tool, status in results.items():
        symbol = "✓" if status else "✗"
        print(f"{symbol} {tool}: {'PASS' if status else 'FAIL'}")
    
    if all(results.values()):
        print("\n✓ All tests passed! Ready to start NTO analysis.")
        return 0
    else:
        print("\n✗ Some tests failed. Please fix issues before proceeding.")
        return 1

if __name__ == '__main__':
    sys.exit(main())
