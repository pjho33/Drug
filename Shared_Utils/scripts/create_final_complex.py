#!/usr/bin/env python3
"""
Create GLUT1-SDG Complex for CHARMM-GUI submission
"""

import os
import shutil

print("=" * 80)
print("Creating GLUT1-SDG Complex")
print("=" * 80)

# File paths
glut1_pdb = '/home/pjho3/바탕화면/GLUT1 only.pdb/charmm-gui-6768295464/step1_pdbreader.pdb'
sdg_pdb = '/home/pjho3/바탕화면/Ligand - 2PEG leg/ligandrm.pdb'
sdg_sdf = '/home/pjho3/바탕화면/Ligand - 2PEG leg/SDG.sdf'
sdg_rtf = '/home/pjho3/바탕화면/Ligand - 2PEG leg/sdg.rtf'
sdg_prm = '/home/pjho3/바탕화면/Ligand - 2PEG leg/sdg.prm'

output_dir = '/home/pjho3/projects/Drug/final_complex'
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

# Step 3: Create complex (GLUT1 + SDG)
print("\n" + "=" * 80)
print("Step 3: Creating GLUT1-SDG complex")
print("=" * 80)

# Renumber SDG atoms
last_atom_num = len(glut1_atoms)
sdg_renumbered = []

for line in sdg_atoms:
    last_atom_num += 1
    # Keep residue name as SDG
    new_line = line[:6] + f"{last_atom_num:>5}" + line[11:]
    sdg_renumbered.append(new_line)

# Combine GLUT1 + SDG
complex_pdb = f'{output_dir}/glut1_sdg_complex.pdb'
with open(complex_pdb, 'w') as f:
    # Write header
    f.write("REMARK   GLUT1-SDG Complex for CHARMM-GUI Membrane Builder\n")
    f.write("REMARK   GLUT1: Glycosylated structure\n")
    f.write("REMARK   SDG: Short-Tripod with 2PEG legs\n")
    
    # Write GLUT1
    for line in glut1_atoms:
        f.write(line)
    
    # Write SDG
    for line in sdg_renumbered:
        f.write(line)
    
    f.write("END\n")

print(f"  ✅ Saved: {complex_pdb}")
print(f"  Total atoms: {len(glut1_atoms) + len(sdg_atoms)}")

# Step 4: Copy SDF file
print("\n" + "=" * 80)
print("Step 4: Copying SDG SDF file")
print("=" * 80)

sdg_sdf_output = f'{output_dir}/sdg.sdf'
shutil.copy(sdg_sdf, sdg_sdf_output)
print(f"  ✅ Copied: {sdg_sdf_output}")

# Step 5: Copy parameter files
print("\n" + "=" * 80)
print("Step 5: Copying parameter files")
print("=" * 80)

shutil.copy(sdg_rtf, f'{output_dir}/sdg.rtf')
shutil.copy(sdg_prm, f'{output_dir}/sdg.prm')
print(f"  ✅ SDG parameters copied")

# Step 6: Verify files
print("\n" + "=" * 80)
print("Step 6: Verification")
print("=" * 80)

print(f"\n  Files in {output_dir}:")
for filename in os.listdir(output_dir):
    filepath = os.path.join(output_dir, filename)
    size = os.path.getsize(filepath)
    print(f"    - {filename}: {size/1024:.1f} KB")

print("\n" + "=" * 80)
print("✅ Complex Generation Complete!")
print("=" * 80)
print(f"\nOutput directory: {output_dir}")
print("\nFiles for CHARMM-GUI submission:")
print("  1. glut1_sdg_complex.pdb - Upload as PDB")
print("  2. sdg.sdf - Upload as Ligand")
print("\nNext steps:")
print("  1. Go to CHARMM-GUI Membrane Builder")
print("  2. Upload glut1_sdg_complex.pdb")
print("  3. Upload sdg.sdf as ligand")
print("  4. Select OpenMM output option ✅")
print("  5. Configure membrane (POPC/POPE)")
print("  6. Submit and download results")
print("=" * 80)
