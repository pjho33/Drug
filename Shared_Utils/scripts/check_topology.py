#!/usr/bin/env python3
"""
Check topology residue information using parmed
"""

import parmed as pmd
from pathlib import Path
from collections import Counter

PARM7_FILE = "/home/pjho3/projects/Drug/final_complex/GLUT1SDGComplex260110/amber/step5_input.parm7"

print("=" * 80)
print("Topology Residue Information")
print("=" * 80)
print()

# Load topology
print(f"Loading: {PARM7_FILE}")
parm = pmd.load_file(PARM7_FILE)
print(f"✅ Loaded successfully")
print()

# Get residue information
print(f"Total atoms: {len(parm.atoms)}")
print(f"Total residues: {len(parm.residues)}")
print()

# Count residue types
res_types = Counter([res.name for res in parm.residues])

print("=" * 80)
print("Residue Type Summary")
print("=" * 80)
print()

# Separate protein, ligand, water, ions
protein_residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
                    'HIE', 'HID', 'HIP', 'CYX', 'NVAL', 'CVAL']
water_residues = ['WAT', 'HOH', 'TIP3', 'TIP4', 'TIP5', 'SPC']
ion_residues = ['NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'SOD', 'CLA', 'POT']

protein_res = {k: v for k, v in res_types.items() if k in protein_residues}
water_res = {k: v for k, v in res_types.items() if k in water_residues}
ion_res = {k: v for k, v in res_types.items() if k in ion_residues}
other_res = {k: v for k, v in res_types.items() if k not in protein_residues + water_residues + ion_residues}

print("Protein residues:")
for res, count in sorted(protein_res.items()):
    print(f"  {res}: {count}")
print(f"  Total protein residues: {sum(protein_res.values())}")
print()

print("Water residues:")
for res, count in sorted(water_res.items()):
    print(f"  {res}: {count}")
print()

print("Ion residues:")
for res, count in sorted(ion_res.items()):
    print(f"  {res}: {count}")
print()

print("⭐ OTHER residues (likely LIGAND):")
for res, count in sorted(other_res.items()):
    print(f"  {res}: {count}")
print()

# Find ligand residue numbers
print("=" * 80)
print("Ligand Residue Details")
print("=" * 80)
print()

for res_name in other_res.keys():
    ligand_residues = [i+1 for i, res in enumerate(parm.residues) if res.name == res_name]
    print(f"Residue name: {res_name}")
    print(f"  Residue numbers (1-indexed): {ligand_residues}")
    print(f"  Total residues: {len(ligand_residues)}")
    
    # Get atom count for first instance
    if ligand_residues:
        first_res_idx = ligand_residues[0] - 1
        atoms_in_res = len(parm.residues[first_res_idx].atoms)
        print(f"  Atoms per residue: {atoms_in_res}")
    print()

# Find protein residue range
protein_res_nums = [i+1 for i, res in enumerate(parm.residues) if res.name in protein_residues]
if protein_res_nums:
    print("=" * 80)
    print("Protein Residue Range")
    print("=" * 80)
    print(f"First protein residue: {min(protein_res_nums)}")
    print(f"Last protein residue: {max(protein_res_nums)}")
    print(f"Total protein residues: {len(protein_res_nums)}")
    print()

# Suggest masks
print("=" * 80)
print("Suggested MMPBSA Masks")
print("=" * 80)
print()

if other_res:
    ligand_name = list(other_res.keys())[0]
    print(f"Ligand mask (by name):  ':{ligand_name}'")
    print(f"Receptor mask (by name): '!:{ligand_name}'")
    print()
    
    ligand_res_list = [i+1 for i, res in enumerate(parm.residues) if res.name == ligand_name]
    if len(ligand_res_list) == 1:
        print(f"Ligand mask (by number): ':{ligand_res_list[0]}'")
        print(f"Receptor mask (by number): '!:{ligand_res_list[0]}'")
    else:
        res_range = f"{min(ligand_res_list)}-{max(ligand_res_list)}"
        print(f"Ligand mask (by range): ':{res_range}'")
        print(f"Receptor mask (by range): '!:{res_range}'")

print()
print("=" * 80)
print("✅ Analysis Complete")
print("=" * 80)
