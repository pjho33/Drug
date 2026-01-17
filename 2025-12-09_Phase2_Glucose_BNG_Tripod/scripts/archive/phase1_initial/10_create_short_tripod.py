#!/usr/bin/env python
"""
Short-Tripod SMILES Generator for Phase 2 Simulation
=====================================================
Design Specs:
- Model A (Basic): TRIS-(PEG2-L-Glucose)3
- Model B (Arginine): Arg-TRIS-(PEG2-L-Glucose)3 (Cation-π interaction)
- PEG length: 1-2 units (entropy minimization)
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import numpy as np
import os

OUTPUT_DIR = "/home/pjho3/projects/Drug/structures"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("=" * 60)
print("Short-Tripod SMILES Generator")
print("=" * 60)

# L-glucose (beta-L-glucopyranose) for O-glycosidic linkage
# Connected via anomeric carbon (C1)
L_GLUCOSE = "O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"

# PEG units
PEG1 = "COCCO"      # 1 ethylene glycol unit
PEG2 = "COCCOCCO"   # 2 ethylene glycol units

def create_model_a_tris(peg_units=2):
    """
    Model A: TRIS-(PEG-L-Glucose)3
    TRIS = tris(hydroxymethyl)aminomethane core: C(CO)(CO)(CO)N
    """
    peg = PEG2 if peg_units == 2 else PEG1
    
    # TRIS core with 3 arms: N-C(CH2-O-PEG-Glucose)3
    # Each arm: -CH2-O-PEG-O-Glucose
    arm = f"CO{peg}{L_GLUCOSE}"
    
    # TRIS: NC(arm)(arm)arm
    smiles = f"NC({arm})({arm}){arm}"
    
    return smiles

def create_model_b_arg_tris(peg_units=2):
    """
    Model B: Arg-TRIS-(PEG-L-Glucose)3
    Arginine provides guanidinium (+) for cation-π interaction
    Structure: Arg-NH-C(CH2-O-PEG-Glucose)3
    """
    peg = PEG2 if peg_units == 2 else PEG1
    
    # Arginine: H2N-C(=NH)-NH-(CH2)3-CH(NH2)-COOH
    # Simplified: guanidinium-propyl chain connected to TRIS nitrogen
    # Arginine side chain: NC(=N)NCCCC
    
    # Full structure: Guanidinium-CH2CH2CH2-NH-TRIS(arms)
    arm = f"CO{peg}{L_GLUCOSE}"
    
    # Arginine-TRIS: [NH2+]=C(N)NCCCNC(arm)(arm)arm
    # Guanidinium is protonated at physiological pH
    smiles = f"NC(=[NH2+])NCCCNC({arm})({arm}){arm}"
    
    return smiles

def create_model_b_simple(peg_units=2):
    """
    Model B simplified: Guanidinium directly on TRIS nitrogen
    """
    peg = PEG2 if peg_units == 2 else PEG1
    arm = f"CO{peg}{L_GLUCOSE}"
    
    # Guanidinium-N-TRIS: NC(=[NH2+])NC(arm)(arm)arm
    smiles = f"NC(=[NH2+])NC({arm})({arm}){arm}"
    
    return smiles

def validate_and_save(smiles, name, output_dir):
    """Validate SMILES and generate 3D structure"""
    print(f"\n{'='*50}")
    print(f"Processing: {name}")
    print(f"{'='*50}")
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"❌ Invalid SMILES for {name}")
        return None
    
    print(f"✅ SMILES valid")
    print(f"   Atoms: {mol.GetNumAtoms()}")
    print(f"   MW: {Descriptors.MolWt(mol):.1f} Da")
    print(f"   Rings: {mol.GetRingInfo().NumRings()}")
    
    # Check for charged groups
    formal_charge = Chem.GetFormalCharge(mol)
    print(f"   Formal charge: {formal_charge:+d}")
    
    # Save SMILES
    smi_path = os.path.join(output_dir, f"{name}.smi")
    with open(smi_path, "w") as f:
        f.write(smiles)
    print(f"   Saved: {smi_path}")
    
    # Generate 3D structure
    print(f"   Generating 3D structure...")
    mol_h = Chem.AddHs(mol)
    
    params = AllChem.ETKDGv3()
    params.maxIterations = 5000
    params.randomSeed = 42
    params.useRandomCoords = True
    
    result = AllChem.EmbedMolecule(mol_h, params)
    if result != 0:
        print(f"   ⚠️ ETKDGv3 failed, trying with random coords...")
        params.useRandomCoords = True
        params.maxIterations = 10000
        result = AllChem.EmbedMolecule(mol_h, params)
    
    if result == 0:
        # Optimize
        try:
            AllChem.MMFFOptimizeMolecule(mol_h, maxIters=2000)
            print(f"   ✅ MMFF optimization done")
        except:
            AllChem.UFFOptimizeMolecule(mol_h, maxIters=2000)
            print(f"   ✅ UFF optimization done")
        
        # Measure size
        conf = mol_h.GetConformer()
        coords = np.array([conf.GetAtomPosition(i) for i in range(mol_h.GetNumAtoms())])
        center = coords.mean(axis=0)
        distances = np.linalg.norm(coords - center, axis=1)
        max_dist = distances.max()
        print(f"   📏 Max radius from center: {max_dist:.1f} Å")
        
        # Save SDF
        sdf_path = os.path.join(output_dir, f"{name}.sdf")
        writer = Chem.SDWriter(sdf_path)
        writer.write(mol_h)
        writer.close()
        print(f"   ✅ Saved: {sdf_path}")
        
        return mol_h
    else:
        print(f"   ❌ 3D embedding failed")
        return None

# Generate structures
print("\n" + "=" * 60)
print("Generating Short-Tripod Structures (PEG2)")
print("=" * 60)

# Model A: Basic TRIS
model_a_smiles = create_model_a_tris(peg_units=2)
print(f"\nModel A SMILES:\n{model_a_smiles}")
mol_a = validate_and_save(model_a_smiles, "model_a_tris_peg2", OUTPUT_DIR)

# Model B: Arginine-TRIS (full arginine chain)
model_b_smiles = create_model_b_arg_tris(peg_units=2)
print(f"\nModel B SMILES:\n{model_b_smiles}")
mol_b = validate_and_save(model_b_smiles, "model_b_arg_tris_peg2", OUTPUT_DIR)

# Model B simplified
model_b_simple_smiles = create_model_b_simple(peg_units=2)
print(f"\nModel B (simple) SMILES:\n{model_b_simple_smiles}")
mol_b_simple = validate_and_save(model_b_simple_smiles, "model_b_guanidinium_tris_peg2", OUTPUT_DIR)

# Also create PEG1 versions for comparison
print("\n" + "=" * 60)
print("Generating Ultra-Short-Tripod Structures (PEG1)")
print("=" * 60)

model_a_peg1 = create_model_a_tris(peg_units=1)
validate_and_save(model_a_peg1, "model_a_tris_peg1", OUTPUT_DIR)

model_b_peg1 = create_model_b_arg_tris(peg_units=1)
validate_and_save(model_b_peg1, "model_b_arg_tris_peg1", OUTPUT_DIR)

print("\n" + "=" * 60)
print("Summary")
print("=" * 60)
print(f"Output directory: {OUTPUT_DIR}")
print("Generated structures:")
print("  - model_a_tris_peg2.sdf (Basic TRIS)")
print("  - model_b_arg_tris_peg2.sdf (Arginine-TRIS, +1 charge)")
print("  - model_b_guanidinium_tris_peg2.sdf (Simplified)")
print("  - model_a_tris_peg1.sdf (Ultra-short)")
print("  - model_b_arg_tris_peg1.sdf (Ultra-short + Arg)")
