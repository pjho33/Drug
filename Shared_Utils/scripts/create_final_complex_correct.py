#!/usr/bin/env python3
"""
Create GLUT1-SDG complex with:
- GLUT1: Desktop clean protein
- SDG structure: Desktop ligand
- SDG position: Phase 2 optimized location
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def format_pdb_line(record, serial, atom_name, res_name, chain_id, res_num, x, y, z, occupancy=1.00, temp=0.00, segment='', element=''):
    """Format proper PDB line"""
    if len(atom_name) == 4:
        atom_field = atom_name
    else:
        atom_field = f" {atom_name:<3}"
    
    line = f"{record:<6}{serial:>5} {atom_field} {res_name:>3} {chain_id:1}{res_num:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}{occupancy:>6.2f}{temp:>6.2f}      {segment:<4}{element:>2}\n"
    return line

print("=" * 80)
print("Creating GLUT1-SDG Complex")
print("GLUT1: Desktop | SDG Structure: Desktop | SDG Position: Phase 2")
print("=" * 80)

# Input files
glut1_file = '/home/pjho3/바탕화면/GLUT1 only.pdb/charmm-gui-6768295464/4pyp_proa.pdb'
desktop_sdf = '/home/pjho3/바탕화면/Ligand - 2PEG leg/SDG.sdf'
phase2_complex = '/home/pjho3/projects/Drug/scripts/Achive/glut1_tripod_complex.pdb'

# Output files
output_pdb = '/home/pjho3/projects/Drug/final_complex/glut1_sdg_complex_optimized_fixed.pdb'
output_sdf = '/home/pjho3/projects/Drug/final_complex/sdg_optimized.sdf'

# ============================================================================
# Step 1: Load desktop SDG structure
# ============================================================================
print("\nStep 1: Loading desktop SDG structure")
mol = Chem.SDMolSupplier(desktop_sdf, removeHs=False)[0]
if mol is None:
    print("❌ Error: Could not read desktop SDF")
    exit(1)

print(f"  ✅ Desktop SDG: {mol.GetNumAtoms()} atoms")
desktop_conf = mol.GetConformer()

# Get desktop SDG center
desktop_coords = []
for i in range(mol.GetNumAtoms()):
    pos = desktop_conf.GetAtomPosition(i)
    desktop_coords.append([pos.x, pos.y, pos.z])
desktop_coords = np.array(desktop_coords)
desktop_center = np.mean(desktop_coords, axis=0)
print(f"  Desktop center: [{desktop_center[0]:.3f}, {desktop_center[1]:.3f}, {desktop_center[2]:.3f}]")

# ============================================================================
# Step 2: Get Phase 2 optimized position
# ============================================================================
print("\nStep 2: Reading Phase 2 optimized position")
with open(phase2_complex, 'r') as f:
    phase2_sdg = [l for l in f if l.startswith('HETATM') and 'SDG' in l]

phase2_coords = []
for line in phase2_sdg:
    parts = line.split()
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
    except:
        continue

phase2_coords = np.array(phase2_coords)
phase2_center = np.mean(phase2_coords, axis=0)
print(f"  ✅ Phase 2 SDG: {len(phase2_coords)} atoms")
print(f"  Phase 2 center: [{phase2_center[0]:.3f}, {phase2_center[1]:.3f}, {phase2_center[2]:.3f}]")

# ============================================================================
# Step 3: Move desktop SDG to Phase 2 position
# ============================================================================
print("\nStep 3: Moving desktop SDG to Phase 2 position")
translation = phase2_center - desktop_center
print(f"  Translation vector: [{translation[0]:.3f}, {translation[1]:.3f}, {translation[2]:.3f}]")

# Apply translation to molecule
for i in range(mol.GetNumAtoms()):
    pos = desktop_conf.GetAtomPosition(i)
    new_pos = (pos.x + translation[0], pos.y + translation[1], pos.z + translation[2])
    desktop_conf.SetAtomPosition(i, new_pos)

# Get new coordinates
moved_coords = []
for i in range(mol.GetNumAtoms()):
    pos = desktop_conf.GetAtomPosition(i)
    moved_coords.append([pos.x, pos.y, pos.z])
moved_coords = np.array(moved_coords)
moved_center = np.mean(moved_coords, axis=0)
print(f"  ✅ Moved center: [{moved_center[0]:.3f}, {moved_center[1]:.3f}, {moved_center[2]:.3f}]")

# ============================================================================
# Step 4: Create complex PDB
# ============================================================================
print("\nStep 4: Creating complex PDB")

# Read GLUT1
with open(glut1_file, 'r') as f:
    glut1_lines = [l for l in f if l.startswith('ATOM')]

pdb_lines = []
pdb_lines.append("REMARK   GLUT1-SDG Complex\n")
pdb_lines.append("REMARK   GLUT1: Desktop 4pyp_proa.pdb (clean, no glycosylation)\n")
pdb_lines.append("REMARK   SDG structure: Desktop SDG.sdf\n")
pdb_lines.append("REMARK   SDG position: Phase 2 optimized location\n")

atom_count = 0

# Write GLUT1
for line in glut1_lines:
    parts = line.split()
    try:
        serial = int(parts[1])
        atom_name = parts[2]
        res_name = parts[3]
        chain_id = parts[4] if parts[4].isalpha() else 'P'
        res_num = int(parts[5]) if parts[4].isalpha() else int(parts[4])
        coord_idx = 6 if parts[4].isalpha() else 5
        
        x = float(parts[coord_idx])
        y = float(parts[coord_idx + 1])
        z = float(parts[coord_idx + 2])
        
        occupancy = float(parts[coord_idx + 3]) if len(parts) > coord_idx + 3 else 1.00
        temp = float(parts[coord_idx + 4]) if len(parts) > coord_idx + 4 else 0.00
        element = atom_name[0]
        
        atom_count += 1
        formatted = format_pdb_line('ATOM', atom_count, atom_name, res_name, chain_id, 
                                   res_num, x, y, z, occupancy, temp, '', element)
        pdb_lines.append(formatted)
    except:
        continue

glut1_count = atom_count
print(f"  GLUT1: {glut1_count} atoms")

# Write SDG with moved coordinates
for i in range(mol.GetNumAtoms()):
    atom = mol.GetAtomWithIdx(i)
    pos = desktop_conf.GetAtomPosition(i)
    atom_name = atom.GetSymbol() + str(i+1)
    element = atom.GetSymbol()
    
    atom_count += 1
    formatted = format_pdb_line('HETATM', atom_count, atom_name, 'SDG', 'L',
                               1, pos.x, pos.y, pos.z, 1.00, 0.00, '', element)
    pdb_lines.append(formatted)

sdg_count = mol.GetNumAtoms()
print(f"  SDG: {sdg_count} atoms")

pdb_lines.append("END\n")

# Write PDB
with open(output_pdb, 'w') as f:
    f.writelines(pdb_lines)
print(f"  ✅ PDB saved: {output_pdb}")

# ============================================================================
# Step 5: Save moved SDG as SDF
# ============================================================================
print("\nStep 5: Saving SDG SDF with Phase 2 position")
writer = Chem.SDWriter(output_sdf)
writer.write(mol)
writer.close()
print(f"  ✅ SDF saved: {output_sdf}")

# ============================================================================
# Verification
# ============================================================================
print("\n" + "=" * 80)
print("Verification")
print("=" * 80)

atom_samples = [l for l in pdb_lines if l.startswith('ATOM')][:3]
hetatm_samples = [l for l in pdb_lines if l.startswith('HETATM')][:3]

print("\nPDB - Sample ATOM lines:")
for line in atom_samples:
    print(f"  {line.rstrip()}")

print("\nPDB - Sample HETATM lines:")
for line in hetatm_samples:
    print(f"  {line.rstrip()}")

verify_mol = Chem.SDMolSupplier(output_sdf, removeHs=False)[0]
if verify_mol:
    verify_conf = verify_mol.GetConformer()
    verify_coords = []
    for i in range(verify_mol.GetNumAtoms()):
        pos = verify_conf.GetAtomPosition(i)
        verify_coords.append([pos.x, pos.y, pos.z])
    verify_center = np.mean(verify_coords, axis=0)
    
    print(f"\nSDF - Atoms: {verify_mol.GetNumAtoms()}")
    print(f"SDF - Center: [{verify_center[0]:.3f}, {verify_center[1]:.3f}, {verify_center[2]:.3f}]")
    print(f"SDF - Sample coordinates:")
    for i in range(min(3, verify_mol.GetNumAtoms())):
        pos = verify_conf.GetAtomPosition(i)
        atom = verify_mol.GetAtomWithIdx(i)
        print(f"  Atom {i+1} ({atom.GetSymbol()}): {pos.x:8.3f} {pos.y:8.3f} {pos.z:8.3f}")

print("\n" + "=" * 80)
print("✅ Files created successfully!")
print("=" * 80)
print(f"\nComplex PDB: {output_pdb}")
print(f"  Total: {atom_count} atoms")
print(f"  GLUT1 (Chain P): {glut1_count} atoms")
print(f"  SDG (Chain L): {sdg_count} atoms")
print(f"\nLigand SDF: {output_sdf}")
print("\nReady for CHARMM-GUI submission!")
