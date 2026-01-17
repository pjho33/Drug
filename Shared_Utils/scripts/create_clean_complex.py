#!/usr/bin/env python3
"""
Create GLUT1-SDG complex using:
- Clean GLUT1 from desktop (no glycosylation)
- SDG ligand from desktop
- Phase 2 optimized position
"""

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
print("Creating Clean GLUT1-SDG Complex")
print("Using desktop files + Phase 2 optimized position")
print("=" * 80)

# Files
glut1_file = '/home/pjho3/바탕화면/GLUT1 only.pdb/charmm-gui-6768295464/4pyp_proa.pdb'
ligand_file = '/home/pjho3/바탕화면/Ligand - 2PEG leg/ligandrm.pdb'
phase2_complex = '/home/pjho3/projects/Drug/scripts/Achive/glut1_tripod_complex.pdb'
output_file = '/home/pjho3/projects/Drug/final_complex/glut1_sdg_complex_optimized_fixed.pdb'

# Step 1: Read clean GLUT1
print("\nStep 1: Reading clean GLUT1 (protein only)")
with open(glut1_file, 'r') as f:
    glut1_lines = [l for l in f if l.startswith('ATOM')]
print(f"  ✅ GLUT1: {len(glut1_lines)} atoms")

# Step 2: Read Phase 2 SDG coordinates
print("\nStep 2: Reading Phase 2 SDG coordinates")
with open(phase2_complex, 'r') as f:
    phase2_sdg_lines = [l for l in f if l.startswith('HETATM') and 'SDG' in l]
print(f"  ✅ Phase 2 SDG: {len(phase2_sdg_lines)} atoms")

# Step 3: Create complex
print("\nStep 3: Creating complex with proper formatting")
output_lines = []
output_lines.append("REMARK   GLUT1 (clean, no glycosylation) + SDG complex\n")
output_lines.append("REMARK   GLUT1: Desktop 4pyp_proa.pdb\n")
output_lines.append("REMARK   SDG: Phase 2 optimized position\n")

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
        output_lines.append(formatted)
    except:
        continue

print(f"  GLUT1 atoms written: {atom_count}")

# Write SDG with Phase 2 coordinates
hetatm_count = 0
for line in phase2_sdg_lines:
    parts = line.split()
    try:
        atom_name = parts[2]
        res_name = 'SDG'
        chain_id = 'L'
        res_num = 1
        
        # Find coordinate index
        coord_idx = 5
        for i, p in enumerate(parts):
            try:
                if '.' in p and float(p):
                    coord_idx = i
                    break
            except:
                continue
        
        x = float(parts[coord_idx])
        y = float(parts[coord_idx + 1])
        z = float(parts[coord_idx + 2])
        
        element = atom_name[0]
        
        atom_count += 1
        hetatm_count += 1
        formatted = format_pdb_line('HETATM', atom_count, atom_name, res_name, chain_id,
                                   res_num, x, y, z, 1.00, 0.00, '', element)
        output_lines.append(formatted)
    except Exception as e:
        print(f"Warning: {e}")
        continue

print(f"  SDG atoms written: {hetatm_count}")

output_lines.append("END\n")

# Write file
print(f"\nWriting: {output_file}")
with open(output_file, 'w') as f:
    f.writelines(output_lines)

print("\n" + "=" * 80)
print("Verification")
print("=" * 80)

atom_samples = [l for l in output_lines if l.startswith('ATOM')][:3]
hetatm_samples = [l for l in output_lines if l.startswith('HETATM')][:3]

print("\nSample ATOM lines:")
for line in atom_samples:
    print(line.rstrip())

print("\nSample HETATM lines:")
for line in hetatm_samples:
    print(line.rstrip())

print("\n" + "=" * 80)
print("✅ Clean complex created!")
print("=" * 80)
print(f"\nTotal atoms: {atom_count}")
print(f"  Protein: {atom_count - hetatm_count}")
print(f"  SDG: {hetatm_count}")
print(f"\nOutput: {output_file}")
