#!/usr/bin/env python3
"""
Fix PDB column alignment for GLUT1-SDG complex
Ensures proper PDB format with Chain ID at column 22
"""

import sys

def fix_pdb_line(line):
    """Fix PDB ATOM/HETATM line to proper column alignment"""
    if not (line.startswith('ATOM') or line.startswith('HETATM')):
        return line
    
    # Parse the line
    record = line[0:6].strip()
    try:
        atom_num = int(line[6:11].strip())
        atom_name = line[12:16].strip()
        alt_loc = line[16:17].strip()
        res_name = line[17:20].strip()
        chain_id = line[21:22].strip()
        res_num = int(line[22:26].strip())
        icode = line[26:27].strip()
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        occupancy = float(line[54:60].strip()) if line[54:60].strip() else 1.00
        temp = float(line[60:66].strip()) if line[60:66].strip() else 0.00
        segment = line[72:76].strip() if len(line) > 72 else ''
        element = line[76:78].strip() if len(line) > 76 else ''
    except (ValueError, IndexError) as e:
        print(f"Warning: Could not parse line: {line.strip()}", file=sys.stderr)
        return line
    
    # Reconstruct with proper alignment
    # PDB format specification:
    # ATOM/HETATM (1-6), serial (7-11), atom name (13-16), altLoc (17),
    # resName (18-20), chainID (22), resSeq (23-26), iCode (27),
    # x (31-38), y (39-46), z (47-54), occupancy (55-60), tempFactor (61-66),
    # segID (73-76), element (77-78)
    
    new_line = f"{record:<6}{atom_num:>5} {atom_name:<4}{alt_loc:1}{res_name:>3} {chain_id:1}{res_num:>4}{icode:1}   {x:>8.3f}{y:>8.3f}{z:>8.3f}{occupancy:>6.2f}{temp:>6.2f}      {segment:<4}{element:>2}\n"
    
    return new_line

def fix_pdb_file(input_file, output_file):
    """Fix entire PDB file"""
    print(f"Reading: {input_file}")
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    fixed_lines = []
    atom_count = 0
    hetatm_count = 0
    
    for line in lines:
        if line.startswith('ATOM'):
            atom_count += 1
            fixed_lines.append(fix_pdb_line(line))
        elif line.startswith('HETATM'):
            hetatm_count += 1
            fixed_lines.append(fix_pdb_line(line))
        else:
            fixed_lines.append(line)
    
    print(f"  ATOM records: {atom_count}")
    print(f"  HETATM records: {hetatm_count}")
    print(f"Writing: {output_file}")
    
    with open(output_file, 'w') as f:
        f.writelines(fixed_lines)
    
    print("✅ PDB file fixed!")

if __name__ == '__main__':
    input_pdb = '/home/pjho3/projects/Drug/scripts/Achive/glut1_tripod_complex.pdb'
    output_pdb = '/home/pjho3/projects/Drug/final_complex/glut1_sdg_complex_optimized_fixed.pdb'
    
    print("=" * 80)
    print("Fixing PDB Column Alignment")
    print("=" * 80)
    
    fix_pdb_file(input_pdb, output_pdb)
    
    print("\n" + "=" * 80)
    print("Verification")
    print("=" * 80)
    
    # Show sample lines
    with open(output_pdb, 'r') as f:
        lines = f.readlines()
        atom_lines = [l for l in lines if l.startswith('ATOM')][:5]
        hetatm_lines = [l for l in lines if l.startswith('HETATM')][:5]
    
    print("\nSample ATOM lines:")
    for line in atom_lines:
        print(line.rstrip())
    
    print("\nSample HETATM lines:")
    for line in hetatm_lines:
        print(line.rstrip())
    
    print("\n" + "=" * 80)
    print("✅ Complete!")
    print("=" * 80)
    print(f"\nFixed file: {output_pdb}")
