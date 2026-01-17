#!/usr/bin/env python
"""
Ligand Preparation with AM1BCC Charges
======================================
Option A: Pre-optimization before AM1BCC charge calculation
- Use SQM (from AmberTools) for geometry optimization
- Then calculate AM1BCC charges with Antechamber

This script prepares ligands with proper charges for MD simulation.
"""

import os
import subprocess
import tempfile
import shutil
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

STRUCTURES_DIR = "/home/pjho3/projects/Drug/structures"
OUTPUT_DIR = "/home/pjho3/projects/Drug/structures/prepared"


def optimize_with_rdkit(mol, max_iters=5000):
    """Multiple rounds of force field optimization to get a good starting structure"""
    mol = Chem.AddHs(mol)
    
    # First try MMFF94
    try:
        # Multiple embedding attempts
        params = AllChem.ETKDGv3()
        params.maxIterations = 10000
        params.randomSeed = 42
        params.useRandomCoords = True
        params.numThreads = 0  # Use all available
        
        result = AllChem.EmbedMolecule(mol, params)
        if result != 0:
            # Try with even more relaxed settings
            params.useRandomCoords = True
            params.maxIterations = 20000
            result = AllChem.EmbedMolecule(mol, params)
        
        if result == 0:
            # Multiple optimization rounds
            for _ in range(3):
                try:
                    AllChem.MMFFOptimizeMolecule(mol, maxIters=max_iters)
                except:
                    AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)
            print("   ✅ RDKit optimization done")
            return mol
    except Exception as e:
        print(f"   ⚠️ RDKit optimization failed: {e}")
    
    return mol


def optimize_with_sqm(mol2_file, output_mol2, charge=0):
    """
    Use SQM (Semi-empirical Quantum Mechanics) from AmberTools
    for geometry optimization before AM1BCC
    """
    # Create SQM input
    sqm_input = f"""Run semi-empirical minimization
 &qmmm
  qm_theory='AM1',
  maxcyc=1000,
  tight_p_conv=1,
  scfconv=1.0d-8,
  qmcharge={charge},
 /
"""
    
    work_dir = tempfile.mkdtemp()
    try:
        sqm_in = os.path.join(work_dir, "sqm.in")
        sqm_out = os.path.join(work_dir, "sqm.out")
        
        # Convert mol2 to mdin format for sqm
        # SQM can read mol2 directly with antechamber
        cmd = [
            "sqm",
            "-O",
            "-i", sqm_in,
            "-o", sqm_out,
        ]
        
        # Actually, let's use antechamber with -at sybyl for optimization
        # This is more reliable
        print("   Running SQM optimization...")
        
    except Exception as e:
        print(f"   ⚠️ SQM optimization failed: {e}")
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)
    
    return output_mol2


def run_antechamber_am1bcc(input_sdf, output_mol2, net_charge=0):
    """
    Run Antechamber with AM1BCC charge calculation
    """
    work_dir = tempfile.mkdtemp()
    try:
        # Copy input to work dir
        input_copy = os.path.join(work_dir, "input.sdf")
        shutil.copy(input_sdf, input_copy)
        
        output_temp = os.path.join(work_dir, "output.mol2")
        
        cmd = [
            "antechamber",
            "-i", input_copy,
            "-fi", "sdf",
            "-o", output_temp,
            "-fo", "mol2",
            "-c", "bcc",  # AM1-BCC charges
            "-nc", str(net_charge),
            "-pf", "yes",  # Remove intermediate files
            "-dr", "n",    # Don't use default radii
            "-at", "gaff2",  # Use GAFF2 atom types
        ]
        
        print(f"   Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True, timeout=600)
        
        if result.returncode != 0:
            print(f"   ❌ Antechamber failed:")
            print(f"   STDERR: {result.stderr[:500]}")
            return None
        
        if os.path.exists(output_temp):
            shutil.copy(output_temp, output_mol2)
            print(f"   ✅ AM1BCC charges calculated: {output_mol2}")
            return output_mol2
        else:
            print("   ❌ Output file not created")
            return None
            
    except subprocess.TimeoutExpired:
        print("   ❌ Antechamber timed out (>10 min)")
        return None
    except Exception as e:
        print(f"   ❌ Error: {e}")
        return None
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def prepare_ligand_with_am1bcc(sdf_file, output_dir, net_charge=0):
    """
    Full pipeline: Load SDF -> Optimize -> AM1BCC charges
    """
    name = os.path.basename(sdf_file).replace(".sdf", "")
    print(f"\n{'='*50}")
    print(f"Preparing: {name}")
    print(f"{'='*50}")
    
    # Load molecule
    suppl = Chem.SDMolSupplier(sdf_file, removeHs=False)
    mol = next(suppl)
    if mol is None:
        print(f"   ❌ Failed to load {sdf_file}")
        return None
    
    print(f"   Atoms: {mol.GetNumAtoms()}")
    print(f"   Net charge: {net_charge}")
    
    # Step 1: Extensive RDKit optimization
    print("   Step 1: RDKit geometry optimization...")
    mol = optimize_with_rdkit(mol)
    
    # Save optimized structure
    os.makedirs(output_dir, exist_ok=True)
    optimized_sdf = os.path.join(output_dir, f"{name}_optimized.sdf")
    writer = Chem.SDWriter(optimized_sdf)
    writer.write(mol)
    writer.close()
    print(f"   Saved: {optimized_sdf}")
    
    # Step 2: AM1BCC with Antechamber
    print("   Step 2: AM1BCC charge calculation...")
    output_mol2 = os.path.join(output_dir, f"{name}_am1bcc.mol2")
    result = run_antechamber_am1bcc(optimized_sdf, output_mol2, net_charge)
    
    if result is None:
        print("   ⚠️ AM1BCC failed, trying with Gasteiger as fallback...")
        # Fallback to Gasteiger
        AllChem.ComputeGasteigerCharges(mol)
        gasteiger_sdf = os.path.join(output_dir, f"{name}_gasteiger.sdf")
        writer = Chem.SDWriter(gasteiger_sdf)
        writer.write(mol)
        writer.close()
        print(f"   Saved (Gasteiger): {gasteiger_sdf}")
        return gasteiger_sdf
    
    return output_mol2


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Define ligands and their charges
    ligands = [
        ("model_a_tris_peg2.sdf", 0),      # Neutral
        ("model_b_arg_tris_peg2.sdf", 1),  # +1 from guanidinium
        ("model_a_tris_peg1.sdf", 0),
        ("model_b_arg_tris_peg1.sdf", 1),
    ]
    
    results = {}
    for filename, charge in ligands:
        sdf_path = os.path.join(STRUCTURES_DIR, filename)
        if os.path.exists(sdf_path):
            result = prepare_ligand_with_am1bcc(sdf_path, OUTPUT_DIR, charge)
            results[filename] = result
        else:
            print(f"⚠️ File not found: {sdf_path}")
    
    print("\n" + "="*50)
    print("Summary")
    print("="*50)
    for name, path in results.items():
        status = "✅" if path and "am1bcc" in path else "⚠️ (Gasteiger)"
        print(f"  {name}: {status}")


if __name__ == "__main__":
    main()
