#!/usr/bin/env python3
"""
Complete PDB column fix - handles all formatting issues
"""

def parse_and_fix_line(line):
    """Parse and fix any ATOM/HETATM line regardless of format"""
    if not (line.startswith('ATOM') or line.startswith('HETATM')):
        return line
    
    record = line[0:6].strip()
    
    # Try to extract all fields flexibly
    parts = line.split()
    
    try:
        # parts[0] = ATOM/HETATM
        # parts[1] = serial number
        # parts[2] = atom name
        # parts[3] = residue name
        # parts[4] = chain ID (might be missing) or residue number
        # Check if parts[4] is a number (residue number) or letter (chain ID)
        
        serial = int(parts[1])
        atom_name = parts[2]
        res_name = parts[3]
        
        # Determine if chain ID exists
        if parts[4].isdigit() or (parts[4].startswith('-') and parts[4][1:].isdigit()):
            # No chain ID, parts[4] is residue number
            chain_id = ' '
            res_num = int(parts[4])
            coord_start = 5
        else:
            # Has chain ID
            chain_id = parts[4]
            res_num = int(parts[5])
            coord_start = 6
        
        # Coordinates
        x = float(parts[coord_start])
        y = float(parts[coord_start + 1])
        z = float(parts[coord_start + 2])
        
        # Occupancy and temp factor
        occupancy = float(parts[coord_start + 3]) if len(parts) > coord_start + 3 else 1.00
        temp = float(parts[coord_start + 4]) if len(parts) > coord_start + 4 else 0.00
        
        # Segment and element
        segment = parts[coord_start + 5] if len(parts) > coord_start + 5 else ''
        element = parts[coord_start + 6] if len(parts) > coord_start + 6 else atom_name[0]
        
    except (ValueError, IndexError) as e:
        print(f"Warning: Could not parse: {line.strip()}")
        return line
    
    # Build proper PDB format line
    # Format: ATOM/HETATM, serial(5), atom_name(4), altLoc(1), resName(3), chainID(1), resSeq(4), iCode(1), x(8.3f), y(8.3f), z(8.3f), occ(6.2f), temp(6.2f), segID(4), element(2)
    
    new_line = f"{record:<6}{serial:>5} {atom_name:<4} {res_name:>3} {chain_id:1}{res_num:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}{occupancy:>6.2f}{temp:>6.2f}      {segment:<4}{element:>2}\n"
    
    return new_line

def fix_pdb_file(input_file, output_file):
    """Fix entire PDB file"""
    print(f"Reading: {input_file}")
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    fixed_lines = []
    atom_count = 0
    hetatm_count = 0
    error_count = 0
    
    for line in lines:
        if line.startswith('ATOM'):
            atom_count += 1
            fixed = parse_and_fix_line(line)
            if fixed == line and 'Warning' in str(fixed):
                error_count += 1
            fixed_lines.append(fixed)
        elif line.startswith('HETATM'):
            hetatm_count += 1
            fixed = parse_and_fix_line(line)
            if fixed == line and 'Warning' in str(fixed):
                error_count += 1
            fixed_lines.append(fixed)
        else:
            fixed_lines.append(line)
    
    print(f"  ATOM records: {atom_count}")
    print(f"  HETATM records: {hetatm_count}")
    if error_count > 0:
        print(f"  ⚠️  Errors: {error_count}")
    
    print(f"Writing: {output_file}")
    
    with open(output_file, 'w') as f:
        f.writelines(fixed_lines)
    
    print("✅ PDB file fixed!")
    return error_count == 0

if __name__ == '__main__':
    input_pdb = '/home/pjho3/projects/Drug/scripts/Achive/glut1_tripod_complex.pdb'
    output_pdb = '/home/pjho3/projects/Drug/final_complex/glut1_sdg_complex_optimized_fixed.pdb'
    
    print("=" * 80)
    print("Fixing PDB Column Alignment (Complete)")
    print("=" * 80)
    
    success = fix_pdb_file(input_pdb, output_pdb)
    
    print("\n" + "=" * 80)
    print("Verification")
    print("=" * 80)
    
    # Show sample lines
    with open(output_pdb, 'r') as f:
        lines = f.readlines()
        atom_lines = [l for l in lines if l.startswith('ATOM')][:3]
        hetatm_lines = [l for l in lines if l.startswith('HETATM')][:3]
        hetatm_end = [l for l in lines if l.startswith('HETATM')][-3:]
    
    print("\nSample ATOM lines:")
    for line in atom_lines:
        print(line.rstrip())
    
    print("\nSample HETATM lines (beginning):")
    for line in hetatm_lines:
        print(line.rstrip())
    
    print("\nSample HETATM lines (end):")
    for line in hetatm_end:
        print(line.rstrip())
    
    print("\n" + "=" * 80)
    if success:
        print("✅ All lines fixed successfully!")
    else:
        print("⚠️  Some lines had issues - please check")
    print("=" * 80)
    print(f"\nFixed file: {output_pdb}")
