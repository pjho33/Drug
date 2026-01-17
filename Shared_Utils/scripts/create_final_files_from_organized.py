#!/usr/bin/env python3
"""
Create both PDB and SDF files together from desktop files
- GLUT1: Desktop clean protein (no glycosylation)
- SDG: Desktop ligand with Phase 2 optimized position
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
print("Creating GLUT1-SDG Complex Files (PDB + SDF)")
print("=" * 80)

# Input files
glut1_file = '/home/pjho3/바탕화면/GLUT1 only.pdb/charmm-gui-6768295464/4pyp_proa.pdb'
phase2_complex = '/home/pjho3/projects/Drug/scripts/Achive/glut1_tripod_complex.pdb'  # Phase 2 optimized position
ligand_sdf = '/home/pjho3/바탕화면/Ligand - 2PEG leg/SDG.sdf'

# Output files
output_pdb = '/home/pjho3/projects/Drug/final_complex/glut1_sdg_complex_optimized_fixed.pdb'
output_sdf = '/home/pjho3/projects/Drug/final_complex/sdg_optimized.sdf'

# ============================================================================
# Part 1: Read ligand PDB to get coordinates
# ============================================================================
print("\nStep 1: Reading Phase 2 optimized SDG coordinates")
with open(phase2_complex, 'r') as f:
    ligand_lines = [l for l in f if l.startswith('HETATM') and 'SDG' in l]

ligand_coords = []
ligand_atoms = []

for line in ligand_lines:
    parts = line.split()
    try:
        atom_name = parts[2]
        
        # Find coordinates
        coord_idx = 5
        for i, p in enumerate(parts):
            try:
                if '.' in p and abs(float(p)) > 0.01:
                    coord_idx = i
                    break
            except:
                continue
        
        x = float(parts[coord_idx])
        y = float(parts[coord_idx + 1])
        z = float(parts[coord_idx + 2])
        
        ligand_coords.append([x, y, z])
        ligand_atoms.append(atom_name)
    except:
        continue

ligand_coords = np.array(ligand_coords)
print(f"  ✅ Ligand: {len(ligand_coords)} atoms")

# ============================================================================
# Part 2: Create PDB file
# ============================================================================
print("\nStep 2: Creating complex PDB")

# Read GLUT1
with open(glut1_file, 'r') as f:
    glut1_lines = [l for l in f if l.startswith('ATOM')]

pdb_lines = []
pdb_lines.append("REMARK   GLUT1 (clean, no glycosylation) + SDG complex\n")
pdb_lines.append("REMARK   GLUT1: Desktop 4pyp_proa.pdb\n")
pdb_lines.append("REMARK   SDG: Phase 2 optimized position\n")

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

# Write SDG
for i, (coord, atom_name) in enumerate(zip(ligand_coords, ligand_atoms)):
    atom_count += 1
    x, y, z = coord
    element = atom_name[0]
    
    formatted = format_pdb_line('HETATM', atom_count, atom_name, 'SDG', 'L',
                               1, x, y, z, 1.00, 0.00, '', element)
    pdb_lines.append(formatted)

sdg_count = len(ligand_coords)
print(f"  SDG: {sdg_count} atoms")

pdb_lines.append("END\n")

# Write PDB
with open(output_pdb, 'w') as f:
    f.writelines(pdb_lines)

print(f"  ✅ PDB saved: {output_pdb}")

# ============================================================================
# Part 3: Create SDF file
# ============================================================================
print("\nStep 3: Creating SDG SDF with matching coordinates")

# Read original SDF structure
mol = Chem.SDMolSupplier(ligand_sdf, removeHs=False)[0]
if mol is None:
    print("  ❌ Error: Could not read original SDF")
else:
    print(f"  Original SDF: {mol.GetNumAtoms()} atoms")
    
    # Update coordinates
    conf = mol.GetConformer()
    
    # Match coordinates (use first N coordinates if counts differ)
    n_coords = min(mol.GetNumAtoms(), len(ligand_coords))
    
    for i in range(n_coords):
        x, y, z = ligand_coords[i]
        conf.SetAtomPosition(i, (float(x), float(y), float(z)))
    
    # Write SDF
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

# PDB verification
atom_samples = [l for l in pdb_lines if l.startswith('ATOM')][:3]
hetatm_samples = [l for l in pdb_lines if l.startswith('HETATM')][:3]

print("\nPDB - Sample ATOM lines:")
for line in atom_samples:
    print(f"  {line.rstrip()}")

print("\nPDB - Sample HETATM lines:")
for line in hetatm_samples:
    print(f"  {line.rstrip()}")

# SDF verification
verify_mol = Chem.SDMolSupplier(output_sdf, removeHs=False)[0]
if verify_mol:
    verify_conf = verify_mol.GetConformer()
    print(f"\nSDF - Atoms: {verify_mol.GetNumAtoms()}")
    print(f"SDF - Sample coordinates:")
    for i in range(min(3, verify_mol.GetNumAtoms())):
        pos = verify_conf.GetAtomPosition(i)
        atom = verify_mol.GetAtomWithIdx(i)
        print(f"  Atom {i+1} ({atom.GetSymbol()}): {pos.x:8.3f} {pos.y:8.3f} {pos.z:8.3f}")

print("\n" + "=" * 80)
print("✅ Both files created successfully!")
print("=" * 80)
print(f"\nComplex PDB: {output_pdb}")
print(f"  Total: {atom_count} atoms")
print(f"  GLUT1 (Chain P): {glut1_count} atoms")
print(f"  SDG (Chain L): {sdg_count} atoms")
print(f"\nLigand SDF: {output_sdf}")
print("\nReady for CHARMM-GUI submission!")
