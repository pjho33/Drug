#!/usr/bin/env python3
"""
Create SDG SDF file from Phase 2 optimized coordinates
Matches the complex PDB file
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

print("=" * 80)
print("Creating SDG SDF from Phase 2 Optimized Coordinates")
print("=" * 80)

# Read Phase 2 SDG coordinates from the complex
phase2_complex = '/home/pjho3/projects/Drug/scripts/Achive/glut1_tripod_complex.pdb'
original_sdf = '/home/pjho3/바탕화면/Ligand - 2PEG leg/SDG.sdf'
output_sdf = '/home/pjho3/projects/Drug/final_complex/sdg_optimized.sdf'

print(f"\nReading original SDF for structure: {original_sdf}")
mol = Chem.SDMolSupplier(original_sdf, removeHs=False)[0]
if mol is None:
    print("Error: Could not read original SDF")
    exit(1)

print(f"  ✅ Molecule loaded: {mol.GetNumAtoms()} atoms")

print(f"\nReading Phase 2 coordinates: {phase2_complex}")
with open(phase2_complex, 'r') as f:
    sdg_lines = [l for l in f if l.startswith('HETATM') and 'SDG' in l]

print(f"  ✅ Found {len(sdg_lines)} SDG atoms")

# Extract coordinates from Phase 2
phase2_coords = []
atom_names = []

for line in sdg_lines:
    parts = line.split()
    atom_name = parts[2]
    
    # Find coordinate index
    coord_idx = 5
    for i, p in enumerate(parts):
        try:
            if '.' in p and abs(float(p)) > 0.01:
                coord_idx = i
                break
        except:
            continue
    
    try:
        x = float(parts[coord_idx])
        y = float(parts[coord_idx + 1])
        z = float(parts[coord_idx + 2])
        phase2_coords.append([x, y, z])
        atom_names.append(atom_name)
    except:
        continue

phase2_coords = np.array(phase2_coords)
print(f"  ✅ Extracted {len(phase2_coords)} coordinates")

# Update molecule coordinates
if len(phase2_coords) != mol.GetNumAtoms():
    print(f"⚠️  Warning: Coordinate count mismatch")
    print(f"   SDF atoms: {mol.GetNumAtoms()}")
    print(f"   PDB atoms: {len(phase2_coords)}")
    
    # Try to match by atom count
    if mol.GetNumAtoms() < len(phase2_coords):
        # SDF might not have hydrogens
        print("   Using first N coordinates to match SDF atom count")
        phase2_coords = phase2_coords[:mol.GetNumAtoms()]

conf = mol.GetConformer()
for i in range(min(mol.GetNumAtoms(), len(phase2_coords))):
    x, y, z = phase2_coords[i]
    conf.SetAtomPosition(i, (float(x), float(y), float(z)))

print(f"\nWriting optimized SDF: {output_sdf}")
writer = Chem.SDWriter(output_sdf)
writer.write(mol)
writer.close()

print("\n" + "=" * 80)
print("Verification")
print("=" * 80)

# Verify the written file
verify_mol = Chem.SDMolSupplier(output_sdf, removeHs=False)[0]
if verify_mol:
    verify_conf = verify_mol.GetConformer()
    print(f"\n✅ SDF file created successfully")
    print(f"   Atoms: {verify_mol.GetNumAtoms()}")
    
    # Show first 3 coordinates
    print(f"\n   Sample coordinates:")
    for i in range(min(3, verify_mol.GetNumAtoms())):
        pos = verify_conf.GetAtomPosition(i)
        atom = verify_mol.GetAtomWithIdx(i)
        print(f"     Atom {i+1} ({atom.GetSymbol()}): {pos.x:8.3f} {pos.y:8.3f} {pos.z:8.3f}")
else:
    print("❌ Error: Could not verify SDF file")

print("\n" + "=" * 80)
print("✅ SDG SDF file created!")
print("=" * 80)
print(f"\nOutput: {output_sdf}")
