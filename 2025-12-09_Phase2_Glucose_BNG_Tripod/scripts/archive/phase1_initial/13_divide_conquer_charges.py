#!/usr/bin/env python
"""
Divide and Conquer AM1BCC Charge Calculation
=============================================
Split large molecules into fragments, calculate charges separately,
then combine for MD simulation.

Strategy:
1. Calculate charges for small fragments:
   - L-glucose (6 atoms heavy)
   - PEG2 linker
   - TRIS core
   - Arginine/Guanidinium
2. Use these as building blocks with pre-calculated charges
"""

import os
import subprocess
import tempfile
import shutil
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import numpy as np

OUTPUT_DIR = "/home/pjho3/projects/Drug/structures/fragments"


def run_antechamber(input_sdf, output_mol2, net_charge=0, timeout=300):
    """Run Antechamber with AM1BCC"""
    work_dir = tempfile.mkdtemp()
    try:
        input_copy = os.path.join(work_dir, "input.sdf")
        shutil.copy(input_sdf, input_copy)
        output_temp = os.path.join(work_dir, "output.mol2")
        
        cmd = [
            "antechamber",
            "-i", input_copy,
            "-fi", "sdf",
            "-o", output_temp,
            "-fo", "mol2",
            "-c", "bcc",
            "-nc", str(net_charge),
            "-pf", "yes",
            "-dr", "n",
            "-at", "gaff2",
        ]
        
        print(f"   Running antechamber (timeout={timeout}s)...")
        result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True, timeout=timeout)
        
        if result.returncode == 0 and os.path.exists(output_temp):
            shutil.copy(output_temp, output_mol2)
            return True
        else:
            print(f"   ❌ Failed: {result.stderr[:200] if result.stderr else 'Unknown error'}")
            return False
            
    except subprocess.TimeoutExpired:
        print(f"   ❌ Timeout after {timeout}s")
        return False
    except Exception as e:
        print(f"   ❌ Error: {e}")
        return False
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def create_and_optimize_fragment(smiles, name):
    """Create 3D structure and optimize"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"   ❌ Invalid SMILES: {smiles}")
        return None
    
    mol = Chem.AddHs(mol)
    
    params = AllChem.ETKDGv3()
    params.maxIterations = 5000
    params.randomSeed = 42
    
    result = AllChem.EmbedMolecule(mol, params)
    if result != 0:
        params.useRandomCoords = True
        result = AllChem.EmbedMolecule(mol, params)
    
    if result == 0:
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
        except:
            AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
    
    return mol


def prepare_fragment(smiles, name, charge, output_dir):
    """Prepare a single fragment with AM1BCC charges"""
    print(f"\n{'='*50}")
    print(f"Fragment: {name}")
    print(f"SMILES: {smiles}")
    print(f"Charge: {charge}")
    print(f"{'='*50}")
    
    mol = create_and_optimize_fragment(smiles, name)
    if mol is None:
        return None
    
    print(f"   Atoms: {mol.GetNumAtoms()} (heavy: {mol.GetNumHeavyAtoms()})")
    print(f"   MW: {Descriptors.MolWt(mol):.1f}")
    
    # Save SDF
    sdf_path = os.path.join(output_dir, f"{name}.sdf")
    writer = Chem.SDWriter(sdf_path)
    writer.write(mol)
    writer.close()
    
    # Run AM1BCC
    mol2_path = os.path.join(output_dir, f"{name}_am1bcc.mol2")
    success = run_antechamber(sdf_path, mol2_path, charge)
    
    if success:
        print(f"   ✅ AM1BCC done: {mol2_path}")
        return mol2_path
    else:
        print(f"   ⚠️ Using Gasteiger fallback")
        AllChem.ComputeGasteigerCharges(mol)
        gasteiger_path = os.path.join(output_dir, f"{name}_gasteiger.sdf")
        writer = Chem.SDWriter(gasteiger_path)
        writer.write(mol)
        writer.close()
        return gasteiger_path


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Define fragments
    # These are the building blocks of our Short-Tripod
    fragments = [
        # L-glucose (beta form, for glycosidic linkage)
        ("OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", "l_glucose", 0),
        
        # PEG2 with terminal OH (simplified)
        ("OCCOCCOCCO", "peg2_linker", 0),
        
        # TRIS core (tris(hydroxymethyl)aminomethane)
        ("NC(CO)(CO)CO", "tris_core", 0),
        
        # Guanidinium (protonated, +1)
        ("NC(=[NH2+])N", "guanidinium", 1),
        
        # Arginine side chain (simplified, +1)
        ("NCCCCNC(=[NH2+])N", "arg_sidechain", 1),
        
        # Single arm: PEG2-O-L-glucose (neutral)
        ("OCCOCCOCCO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O", "peg2_glucose_arm", 0),
        
        # TRIS with one arm attached (test)
        ("NC(CO)(CO)COCCOCCOCCO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O", "tris_one_arm", 0),
    ]
    
    results = {}
    for smiles, name, charge in fragments:
        result = prepare_fragment(smiles, name, charge, OUTPUT_DIR)
        results[name] = result
    
    # Summary
    print("\n" + "="*60)
    print("Fragment Preparation Summary")
    print("="*60)
    
    success_count = 0
    for name, path in results.items():
        if path:
            method = "AM1BCC" if "am1bcc" in path else "Gasteiger"
            print(f"  ✅ {name}: {method}")
            success_count += 1
        else:
            print(f"  ❌ {name}: FAILED")
    
    print(f"\nSuccess: {success_count}/{len(fragments)}")
    print(f"Output: {OUTPUT_DIR}")
    
    # If fragments work, we can build the full molecule
    if success_count >= 5:
        print("\n" + "="*60)
        print("Next Step: Build full molecule from fragments")
        print("="*60)
        print("Use tleap or similar to combine fragments with their charges")


if __name__ == "__main__":
    main()
