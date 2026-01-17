#!/usr/bin/env python3
"""
Create GLUT1-Ligand Complexes (Experimental and Control)
- Experimental: GLUT1 + SDG (2PEG leg)
- Control: GLUT1 + D-Glucose
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import os

print("=" * 80)
print("Creating GLUT1-Ligand Complexes")
print("=" * 80)

# File paths
glut1_pdb = '/home/pjho3/바탕화면/GLUT1 only.pdb/charmm-gui-6768295464/step1_pdbreader.pdb'
sdg_pdb = '/home/pjho3/바탕화면/Ligand - 2PEG leg/ligandrm.pdb'
sdg_sdf = '/home/pjho3/바탕화면/Ligand - 2PEG leg/SDG.sdf'

output_dir = '/home/pjho3/projects/Drug/complexes'
os.makedirs(output_dir, exist_ok=True)

print(f"\n📁 Output directory: {output_dir}")

# Step 1: Load GLUT1
print("\n" + "=" * 80)
print("Step 1: Loading GLUT1")
print("=" * 80)

with open(glut1_pdb, 'r') as f:
    glut1_lines = f.readlines()

glut1_atoms = [l for l in glut1_lines if l.startswith('ATOM') or l.startswith('HETATM')]
print(f"  ✅ GLUT1 loaded: {len(glut1_atoms)} atoms")

# Step 2: Load SDG ligand
print("\n" + "=" * 80)
print("Step 2: Loading SDG ligand")
print("=" * 80)

with open(sdg_pdb, 'r') as f:
    sdg_lines = f.readlines()

sdg_atoms = [l for l in sdg_lines if l.startswith('ATOM') or l.startswith('HETATM')]
print(f"  ✅ SDG loaded: {len(sdg_atoms)} atoms")

# Step 3: Create experimental complex (GLUT1 + SDG)
print("\n" + "=" * 80)
print("Step 3: Creating experimental complex (GLUT1 + SDG)")
print("=" * 80)

# Renumber SDG atoms
last_atom_num = len(glut1_atoms)
sdg_renumbered = []

for line in sdg_atoms:
    last_atom_num += 1
    # Change residue name to SDG if not already
    new_line = line[:6] + f"{last_atom_num:>5}" + line[11:17] + "SDG " + line[21:]
    sdg_renumbered.append(new_line)

# Combine GLUT1 + SDG
experimental_complex = f'{output_dir}/glut1_sdg_complex.pdb'
with open(experimental_complex, 'w') as f:
    # Write header
    f.write("REMARK   Experimental Complex: GLUT1 + SDG (2PEG leg)\n")
    f.write("REMARK   Generated for Phase 3 MD simulation\n")
    
    # Write GLUT1
    for line in glut1_atoms:
        f.write(line)
    
    # Write SDG
    for line in sdg_renumbered:
        f.write(line)
    
    f.write("END\n")

print(f"  ✅ Saved: {experimental_complex}")
print(f"  Total atoms: {len(glut1_atoms) + len(sdg_atoms)}")

# Step 4: Create control complex (GLUT1 + D-Glucose)
print("\n" + "=" * 80)
print("Step 4: Creating control complex (GLUT1 + D-Glucose)")
print("=" * 80)

# Create D-Glucose from SMILES
glucose_smiles = "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O"  # D-Glucose
print(f"  D-Glucose SMILES: {glucose_smiles}")

mol = Chem.MolFromSmiles(glucose_smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, randomSeed=42)
AllChem.MMFFOptimizeMolecule(mol)

print(f"  ✅ D-Glucose generated: {mol.GetNumAtoms()} atoms")

# Position D-Glucose at the same location as SDG (approximate)
# Extract SDG center of mass
sdg_coords = []
for line in sdg_atoms:
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        sdg_coords.append([x, y, z])
    except:
        continue

sdg_center = np.mean(sdg_coords, axis=0)
print(f"  SDG center: {sdg_center}")

# Move glucose to SDG position
conf = mol.GetConformer()
glucose_coords = conf.GetPositions()
glucose_center = np.mean(glucose_coords, axis=0)
translation = sdg_center - glucose_center

for i in range(mol.GetNumAtoms()):
    pos = conf.GetAtomPosition(i)
    conf.SetAtomPosition(i, (pos.x + translation[0], 
                             pos.y + translation[1], 
                             pos.z + translation[2]))

# Write control complex
control_complex = f'{output_dir}/glut1_glucose_complex.pdb'
with open(control_complex, 'w') as f:
    # Write header
    f.write("REMARK   Control Complex: GLUT1 + D-Glucose\n")
    f.write("REMARK   Generated for Phase 3 MD simulation\n")
    
    # Write GLUT1
    for line in glut1_atoms:
        f.write(line)
    
    # Write D-Glucose
    atom_num = len(glut1_atoms)
    for i, atom in enumerate(mol.GetAtoms()):
        atom_num += 1
        pos = conf.GetAtomPosition(i)
        element = atom.GetSymbol()
        
        pdb_line = f"HETATM{atom_num:>5}  {element:<3} GLC     1    {pos.x:>8.3f}{pos.y:>8.3f}{pos.z:>8.3f}  1.00  0.00           {element:>2}\n"
        f.write(pdb_line)
    
    f.write("END\n")

print(f"  ✅ Saved: {control_complex}")
print(f"  Total atoms: {len(glut1_atoms) + mol.GetNumAtoms()}")

# Step 5: Create/Copy SDF files
print("\n" + "=" * 80)
print("Step 5: Creating SDF files")
print("=" * 80)

# Copy SDG SDF
import shutil
sdg_sdf_output = f'{output_dir}/sdg.sdf'
shutil.copy(sdg_sdf, sdg_sdf_output)
print(f"  ✅ SDG SDF: {sdg_sdf_output}")

# Create D-Glucose SDF
glucose_sdf_output = f'{output_dir}/glucose.sdf'
writer = Chem.SDWriter(glucose_sdf_output)
writer.write(mol)
writer.close()
print(f"  ✅ D-Glucose SDF: {glucose_sdf_output}")

# Step 6: Copy parameter files
print("\n" + "=" * 80)
print("Step 6: Copying parameter files")
print("=" * 80)

sdg_rtf = '/home/pjho3/바탕화면/Ligand - 2PEG leg/sdg.rtf'
sdg_prm = '/home/pjho3/바탕화면/Ligand - 2PEG leg/sdg.prm'

shutil.copy(sdg_rtf, f'{output_dir}/sdg.rtf')
shutil.copy(sdg_prm, f'{output_dir}/sdg.prm')
print(f"  ✅ SDG parameters copied")

print("\n" + "=" * 80)
print("✅ Complex Generation Complete!")
print("=" * 80)
print(f"\nOutput files in: {output_dir}")
print("\nExperimental (GLUT1 + SDG):")
print(f"  - glut1_sdg_complex.pdb")
print(f"  - sdg.sdf")
print(f"  - sdg.rtf, sdg.prm")
print("\nControl (GLUT1 + D-Glucose):")
print(f"  - glut1_glucose_complex.pdb")
print(f"  - glucose.sdf")
print("\n" + "=" * 80)
print("Next steps:")
print("  1. Submit glut1_sdg_complex.pdb + sdg.sdf to CHARMM-GUI")
print("  2. Submit glut1_glucose_complex.pdb + glucose.sdf to CHARMM-GUI")
print("  3. Run MD simulations on both systems")
print("=" * 80)
