#!/usr/bin/env python3
"""
Simple PDB combination: CHARMM-GUI membrane + SDG ligand
"""

print("=" * 80)
print("Combining CHARMM-GUI Membrane + SDG Ligand")
print("=" * 80)

# Read membrane system PDB
membrane_pdb = '/home/pjho3/다운로드/charmm-gui-6750265216membranebuilder/openmm/step5_input.pdb'
sdg_pdb = 'glut1_tripod_complex_sdg_only.pdb'
output_pdb = 'glut1_sdg_membrane_combined.pdb'

print(f"\n1. Reading membrane PDB...")
with open(membrane_pdb, 'r') as f:
    membrane_lines = f.readlines()

membrane_atoms = [l for l in membrane_lines if l.startswith('ATOM') or l.startswith('HETATM')]
other_lines = [l for l in membrane_lines if not (l.startswith('ATOM') or l.startswith('HETATM') or l.startswith('END'))]

print(f"   Membrane atoms: {len(membrane_atoms)}")

print(f"\n2. Reading SDG PDB...")
with open(sdg_pdb, 'r') as f:
    sdg_lines = f.readlines()

sdg_atoms = [l for l in sdg_lines if l.startswith('ATOM') or l.startswith('HETATM')]
print(f"   SDG atoms: {len(sdg_atoms)}")

# Renumber SDG atoms
print(f"\n3. Renumbering SDG atoms...")
last_atom_num = len(membrane_atoms)
renumbered_sdg = []

for line in sdg_atoms:
    last_atom_num += 1
    # Replace atom serial number
    new_line = line[:6] + f"{last_atom_num:>5}" + line[11:]
    renumbered_sdg.append(new_line)

print(f"   SDG atoms renumbered starting from {len(membrane_atoms) + 1}")

# Combine
print(f"\n4. Combining structures...")
with open(output_pdb, 'w') as f:
    # Write header lines
    for line in other_lines:
        f.write(line)
    
    # Write membrane atoms
    for line in membrane_atoms:
        f.write(line)
    
    # Write SDG atoms
    for line in renumbered_sdg:
        f.write(line)
    
    # Write END
    f.write('END\n')

print(f"   ✅ Saved: {output_pdb}")
print(f"   Total atoms: {len(membrane_atoms) + len(sdg_atoms)}")

print("\n" + "=" * 80)
print("✅ PDB Combination Complete!")
print("=" * 80)
print(f"\nOutput: {output_pdb}")
print("\nNext: Run OpenMM with CharmmParameterSet")
print("=" * 80)
