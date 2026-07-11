#!/usr/bin/env python3
"""
Phase 15 Expansion: Orbital Composition Analysis for 17 molecules
Uses Multiwfn to extract atomic contributions to HOMO/LUMO

FIXED VERSION: Implements proper Multiwfn workflow
1. Move calculation folder to /home/tchapet
2. Run Multiwfn in home directory
3. Move folder back to original location
"""

import subprocess
import shutil
from pathlib import Path
import pandas as pd
import json
import re
from datetime import datetime

# All 17 successfully optimized molecules
TARGET_MOLECULES = [10132, 12635, 14380, 17851, 18985, 19531, 20722, 20778, 
                    24255, 24524, 4550, 6678, 6775, 6966, 8371, 8521, 977]

BASE_DIR = Path(__file__).parent.parent.resolve()
INPUT_DIR = BASE_DIR / 'data' / 'optimized_input'
OUTPUT_DIR = BASE_DIR / 'analysis'
WORK_DIR = BASE_DIR / 'data' / 'orbital_expanded'
HOME_DIR = Path('/home/tchapet')
WORK_DIR.mkdir(parents=True, exist_ok=True)
CHECKPOINT_FILE = OUTPUT_DIR / 'orbital_expanded_checkpoint.json'

def load_checkpoint():
    """Load checkpoint to resume interrupted calculations"""
    if CHECKPOINT_FILE.exists():
        with open(CHECKPOINT_FILE) as f:
            return json.load(f)
    return {'completed': [], 'results': []}

def save_checkpoint(data):
    """Save checkpoint after each successful calculation"""
    with open(CHECKPOINT_FILE, 'w') as f:
        json.dump(data, f, indent=2)
    print(f"  [Checkpoint saved: {len(data['completed'])} completed]")

def run_xtb_wavefunction(xyz_file, work_dir):
    """Generate wavefunction with xTB"""
    result = subprocess.run(
        ['xtb', str(xyz_file), '--sp', '--gfn', '2', '--molden'],
        cwd=work_dir,
        capture_output=True,
        text=True,
        timeout=300
    )
    
    molden = work_dir / 'molden.input'
    return molden.exists()

def parse_multiwfn_output(output_text):
    """Parse Multiwfn orbital composition output"""
    lines = output_text.split('\n')
    
    compositions = {}
    current_orbital = None
    in_atom_section = False
    
    for i, line in enumerate(lines):
        # Detect which orbital we're analyzing
        if 'Orbital:' in line:
            if 'Occ:  2.000000' in line or 'Occ:  1.000000' in line:
                # This is likely HOMO (occupied)
                current_orbital = 'HOMO'
                in_atom_section = False
            elif 'Occ:  0.000000' in line:
                # This is LUMO (unoccupied)
                current_orbital = 'LUMO'
                in_atom_section = False
        
        # Detect start of atom composition section
        if 'Composition of each atom:' in line:
            in_atom_section = True
            continue
        
        # Stop at empty line or next section
        if in_atom_section and (line.strip() == '' or 'Orbital delocalization' in line):
            in_atom_section = False
            continue
        
        # Parse atomic contributions (format: " Atom     1(O ) :    28.78472 %")
        if current_orbital and in_atom_section and 'Atom' in line and '%' in line:
            match = re.search(r'Atom\s+\d+\(([A-Z][a-z]?)\s*\)\s*:\s+([\d.]+)\s*%', line)
            if match:
                atom_sym, percent = match.groups()
                key = f"{current_orbital}_{atom_sym}"
                if key not in compositions:
                    compositions[key] = 0.0
                compositions[key] += float(percent)
    
    return compositions

def analyze_orbital_composition(mol_id):
    """
    Analyze orbital composition using Multiwfn with proper home directory workflow
    
    Steps:
    1. Generate wavefunction with xTB in work directory
    2. Move work directory to /home/tchapet
    3. Run Multiwfn in home directory
    4. Move work directory back to original location
    """
    
    print(f"\n[{datetime.now().strftime('%H:%M:%S')}] Processing molecule {mol_id}...")
    
    xyz_file = INPUT_DIR / f"{mol_id}/{mol_id}_optimized.xyz"
    if not xyz_file.exists():
        print(f"✗ {mol_id}: XYZ file not found at {xyz_file}")
        return None
    
    # Create work directory in project
    mol_work = WORK_DIR / str(mol_id)
    mol_work.mkdir(exist_ok=True)
    
    # Generate wavefunction
    print(f"  Generating wavefunction...", end='', flush=True)
    if not run_xtb_wavefunction(xyz_file, mol_work):
        print(" FAILED")
        return None
    print(" done")
    
    # HOME DIRECTORY WORKFLOW STARTS HERE
    home_work = HOME_DIR / f'orbital_temp_{mol_id}'
    
    try:
        # Step 1: Move to home directory
        print(f"  Moving to home directory...", end='', flush=True)
        if home_work.exists():
            shutil.rmtree(home_work)
        shutil.move(str(mol_work), str(home_work))
        print(" done")
        
        molden_file = home_work / 'molden.input'
        if not molden_file.exists():
            print(f"✗ {mol_id}: molden.input not found after move")
            shutil.move(str(home_work), str(mol_work))  # Move back
            return None
        
        # Step 2: Run Multiwfn in home directory
        print(f"  Running Multiwfn in home...", end='', flush=True)
        
        # Multiwfn input commands:
        # 8 = Orbital composition analysis
        # 1 = Mulliken partition analysis
        # h = HOMO
        # l = LUMO
        # 0 = Return to main menu
        # q = Quit
        multiwfn_input = "8\n1\nh\nl\n0\nq\n"
        
        result = subprocess.run(
            ['Multiwfn', str(molden_file)],
            input=multiwfn_input,
            cwd=home_work,
            capture_output=True,
            text=True,
            timeout=120
        )
        
        # Save Multiwfn output log
        log_file = home_work / f'multiwfn_orbital_{mol_id}.log'
        with open(log_file, 'w') as f:
            f.write(result.stdout)
        
        print(" done")
        
        # Step 3: Parse results
        print(f"  Parsing results...", end='', flush=True)
        compositions = parse_multiwfn_output(result.stdout)
        
        if compositions:
            print(" done")
            
            # Extract key contributions
            homo_n = compositions.get('HOMO_N', 0.0)
            lumo_n = compositions.get('LUMO_N', 0.0)
            homo_c = compositions.get('HOMO_C', 0.0)
            lumo_c = compositions.get('LUMO_C', 0.0)
            homo_o = compositions.get('HOMO_O', 0.0)
            lumo_o = compositions.get('LUMO_O', 0.0)
            
            print(f"✓ {mol_id}: HOMO_N={homo_n:.1f}%, LUMO_N={lumo_n:.1f}%")
            
            result_dict = {
                'mol_id': mol_id,
                'HOMO_N_percent': homo_n,
                'LUMO_N_percent': lumo_n,
                'HOMO_C_percent': homo_c,
                'LUMO_C_percent': lumo_c,
                'HOMO_O_percent': homo_o,
                'LUMO_O_percent': lumo_o
            }
            
            # Add all other elements found
            for key, val in compositions.items():
                col_name = f"{key}_percent"
                if col_name not in result_dict:
                    result_dict[col_name] = val
            
            # Step 4: Move back to original location
            print(f"  Moving back to project...", end='', flush=True)
            if mol_work.exists():
                shutil.rmtree(mol_work)
            shutil.move(str(home_work), str(mol_work))
            print(" done")
            
            return result_dict
        else:
            print(" NO DATA")
            # Move back even if parsing failed
            if mol_work.exists():
                shutil.rmtree(mol_work)
            shutil.move(str(home_work), str(mol_work))
            return None
            
    except Exception as e:
        print(f" ERROR: {e}")
        # Ensure we move back even on error
        try:
            if home_work.exists():
                if mol_work.exists():
                    shutil.rmtree(mol_work)
                shutil.move(str(home_work), str(mol_work))
                print(f"  [Folder moved back after error]")
        except Exception as move_error:
            print(f"  [WARNING: Could not move folder back: {move_error}]")
        return None

def main():
    print("="*80)
    print("PHASE 15 EXPANDED: ORBITAL COMPOSITION ANALYSIS (FIXED)")
    print("="*80)
    print(f"Target molecules: {len(TARGET_MOLECULES)}")
    print(f"Method: GFN2-xTB + Multiwfn (with home directory workflow)")
    print(f"Work directory: {WORK_DIR}")
    print(f"Home temp directory: {HOME_DIR}/orbital_temp_*\n")
    
    checkpoint = load_checkpoint()
    completed = set(checkpoint['completed'])
    results = checkpoint['results']
    
    if completed:
        print(f"Resuming: {len(completed)} molecules already completed")
        print(f"Completed: {sorted(completed)}\n")
    
    success_count = 0
    fail_count = 0
    
    for i, mol_id in enumerate(TARGET_MOLECULES, 1):
        print(f"\n[{i}/{len(TARGET_MOLECULES)}]", end=' ')
        
        if mol_id in completed:
            print(f"Skipping {mol_id} (already completed)")
            success_count += 1
            continue
        
        try:
            result = analyze_orbital_composition(mol_id)
            if result:
                results.append(result)
                completed.add(mol_id)
                save_checkpoint({'completed': list(completed), 'results': results})
                success_count += 1
            else:
                fail_count += 1
        except Exception as e:
            print(f"✗ {mol_id}: Unexpected error - {e}")
            fail_count += 1
            continue
    
    # Save final results
    print(f"\n{'='*80}")
    print(f"ANALYSIS COMPLETE")
    print(f"{'='*80}")
    print(f"Success: {success_count}/{len(TARGET_MOLECULES)}")
    print(f"Failed: {fail_count}/{len(TARGET_MOLECULES)}")
    
    if results:
        df = pd.DataFrame(results)
        output_file = OUTPUT_DIR / 'orbital_composition_expanded.csv'
        df.to_csv(output_file, index=False)
        print(f"\n✓ Results saved to {output_file}")
        print(f"✓ Total analyzed: {len(results)} molecules")
        
        # Statistics
        if 'HOMO_N_percent' in df.columns and 'LUMO_N_percent' in df.columns:
            print(f"\n{'='*80}")
            print(f"NITROGEN CONTRIBUTION STATISTICS")
            print(f"{'='*80}")
            print(f"HOMO_N: {df['HOMO_N_percent'].mean():.1f} ± {df['HOMO_N_percent'].std():.1f}% (n={len(df)})")
            print(f"LUMO_N: {df['LUMO_N_percent'].mean():.1f} ± {df['LUMO_N_percent'].std():.1f}% (n={len(df)})")
            print(f"\nMolecules with N-rich HOMO (>50%): {(df['HOMO_N_percent'] > 50).sum()}")
            print(f"Molecules with N-rich LUMO (>50%): {(df['LUMO_N_percent'] > 50).sum()}")
        
        print(f"\n✓ Phase 15 expansion complete!")
        print(f"✓ Sample size increased from 2 → {len(results)} molecules")
        
        # Clean up checkpoint
        if CHECKPOINT_FILE.exists():
            CHECKPOINT_FILE.unlink()
            print(f"✓ Checkpoint file cleaned up")
    else:
        print("\n✗ No results to save")
        print("Check error messages above for details")

if __name__ == '__main__':
    main()
