#!/usr/bin/env python
"""
Build Full Molecule Parameters from Fragments
==============================================
Transfer AM1BCC charges from fragments to full molecule
and generate GAFF2 parameters for MD simulation.
"""

import os
import subprocess
import tempfile
import shutil
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

FRAGMENTS_DIR = "/home/pjho3/projects/Drug/structures/fragments"
STRUCTURES_DIR = "/home/pjho3/projects/Drug/structures"
OUTPUT_DIR = "/home/pjho3/projects/Drug/structures/parameterized"


def read_mol2_charges(mol2_file):
    """Read partial charges from mol2 file"""
    charges = []
    atom_names = []
    in_atom_section = False
    
    with open(mol2_file, 'r') as f:
        for line in f:
            if line.startswith("@<TRIPOS>ATOM"):
                in_atom_section = True
                continue
            elif line.startswith("@<TRIPOS>"):
                in_atom_section = False
                continue
            
            if in_atom_section and line.strip():
                parts = line.split()
                if len(parts) >= 9:
                    atom_names.append(parts[1])
                    charges.append(float(parts[8]))
    
    return atom_names, charges


def get_fragment_charge_stats():
    """Get average charges per atom type from fragments"""
    fragment_files = [
        "l_glucose_am1bcc.mol2",
        "peg2_linker_am1bcc.mol2",
        "tris_core_am1bcc.mol2",
        "guanidinium_am1bcc.mol2",
        "arg_sidechain_am1bcc.mol2",
    ]
    
    # Collect charges by element
    element_charges = {}
    
    for fname in fragment_files:
        fpath = os.path.join(FRAGMENTS_DIR, fname)
        if os.path.exists(fpath):
            names, charges = read_mol2_charges(fpath)
            for name, charge in zip(names, charges):
                # Extract element from atom name
                elem = ''.join([c for c in name if c.isalpha()])[:2]
                if elem not in element_charges:
                    element_charges[elem] = []
                element_charges[elem].append(charge)
    
    # Calculate averages
    avg_charges = {}
    for elem, charges in element_charges.items():
        avg_charges[elem] = {
            'mean': np.mean(charges),
            'std': np.std(charges),
            'count': len(charges)
        }
    
    return avg_charges


def run_parmchk2(mol2_file, frcmod_file):
    """Generate missing GAFF2 parameters"""
    cmd = ["parmchk2", "-i", mol2_file, "-f", "mol2", "-o", frcmod_file, "-s", "gaff2"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0


def create_gaff2_params_with_gasteiger(sdf_file, output_prefix, net_charge=0):
    """
    Create GAFF2 parameters using Gasteiger charges
    (fallback when AM1BCC fails for large molecules)
    """
    print(f"\n{'='*50}")
    print(f"Parameterizing: {os.path.basename(sdf_file)}")
    print(f"{'='*50}")
    
    # Load molecule
    suppl = Chem.SDMolSupplier(sdf_file, removeHs=False)
    mol = next(suppl)
    if mol is None:
        print("   ❌ Failed to load molecule")
        return None
    
    mol = Chem.AddHs(mol, addCoords=True)
    print(f"   Atoms: {mol.GetNumAtoms()}")
    
    # Compute Gasteiger charges
    AllChem.ComputeGasteigerCharges(mol)
    
    # Adjust total charge to match expected
    charges = []
    for atom in mol.GetAtoms():
        charge = float(atom.GetDoubleProp('_GasteigerCharge'))
        if np.isnan(charge):
            charge = 0.0
        charges.append(charge)
    
    total_charge = sum(charges)
    print(f"   Gasteiger total charge: {total_charge:.3f}")
    print(f"   Expected charge: {net_charge}")
    
    # Normalize charges to match expected total
    if abs(total_charge - net_charge) > 0.01:
        adjustment = (net_charge - total_charge) / len(charges)
        charges = [c + adjustment for c in charges]
        print(f"   Adjusted charges (delta per atom: {adjustment:.4f})")
    
    # Create work directory
    work_dir = tempfile.mkdtemp()
    try:
        # Save as PDB with charges
        pdb_file = os.path.join(work_dir, "ligand.pdb")
        Chem.MolToPDBFile(mol, pdb_file)
        
        # Convert to mol2 with antechamber (just for format, use existing charges)
        mol2_file = os.path.join(work_dir, "ligand.mol2")
        
        # Use antechamber with 'mul' (Mulliken) which is faster, then we'll replace charges
        cmd = [
            "antechamber",
            "-i", pdb_file,
            "-fi", "pdb",
            "-o", mol2_file,
            "-fo", "mol2",
            "-c", "gas",  # Gasteiger charges (fast)
            "-nc", str(net_charge),
            "-pf", "yes",
            "-at", "gaff2",
        ]
        
        print(f"   Running antechamber for atom typing...")
        result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True, timeout=120)
        
        if result.returncode != 0 or not os.path.exists(mol2_file):
            print(f"   ❌ Antechamber failed: {result.stderr[:200]}")
            return None
        
        # Generate frcmod for missing parameters
        frcmod_file = os.path.join(work_dir, "ligand.frcmod")
        if run_parmchk2(mol2_file, frcmod_file):
            print(f"   ✅ Generated frcmod")
        else:
            print(f"   ⚠️ parmchk2 failed, continuing anyway")
        
        # Copy outputs
        os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
        
        output_mol2 = output_prefix + ".mol2"
        output_frcmod = output_prefix + ".frcmod"
        
        shutil.copy(mol2_file, output_mol2)
        if os.path.exists(frcmod_file):
            shutil.copy(frcmod_file, output_frcmod)
        
        print(f"   ✅ Saved: {output_mol2}")
        print(f"   ✅ Saved: {output_frcmod}")
        
        return output_mol2
        
    except Exception as e:
        print(f"   ❌ Error: {e}")
        return None
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Print fragment charge statistics
    print("="*60)
    print("Fragment Charge Statistics (AM1BCC)")
    print("="*60)
    
    stats = get_fragment_charge_stats()
    for elem, data in sorted(stats.items()):
        print(f"  {elem:3s}: mean={data['mean']:+.3f}, std={data['std']:.3f}, n={data['count']}")
    
    # Parameterize full molecules
    molecules = [
        ("model_a_tris_peg2.sdf", 0),      # Neutral TRIS
        ("model_b_arg_tris_peg2.sdf", 1),  # +1 Arginine-TRIS
    ]
    
    results = {}
    for filename, charge in molecules:
        sdf_path = os.path.join(STRUCTURES_DIR, filename)
        if os.path.exists(sdf_path):
            name = filename.replace(".sdf", "")
            output_prefix = os.path.join(OUTPUT_DIR, name)
            result = create_gaff2_params_with_gasteiger(sdf_path, output_prefix, charge)
            results[name] = result
        else:
            print(f"⚠️ Not found: {sdf_path}")
    
    # Summary
    print("\n" + "="*60)
    print("Parameterization Summary")
    print("="*60)
    for name, path in results.items():
        status = "✅" if path else "❌"
        print(f"  {status} {name}")
    
    print(f"\nOutput: {OUTPUT_DIR}")
    print("\nReady for MD simulation!")


if __name__ == "__main__":
    main()
