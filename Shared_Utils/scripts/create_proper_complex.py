#!/usr/bin/env python3
"""
Create properly formatted GLUT1-SDG complex with correct PDB columns
Uses Phase 2 optimized positions
"""

def format_pdb_line(record, serial, atom_name, res_name, chain_id, res_num, x, y, z, occupancy=1.00, temp=0.00, segment='', element=''):
    """Format a proper PDB ATOM/HETATM line"""
    # PDB format specification:
    # cols 1-6: record name
    # cols 7-11: atom serial number (right justified)
    # col 12: blank
    # cols 13-16: atom name (left justified for 1-3 char, centered for 4 char)
    # col 17: alternate location indicator
    # cols 18-20: residue name (right justified)
    # col 21: blank
    # col 22: chain identifier
    # cols 23-26: residue sequence number (right justified)
    # col 27: insertion code
    # cols 28-30: blank
    # cols 31-38: X coordinate (8.3f)
    # cols 39-46: Y coordinate (8.3f)
    # cols 47-54: Z coordinate (8.3f)
    # cols 55-60: occupancy (6.2f)
    # cols 61-66: temperature factor (6.2f)
    # cols 67-72: blank
    # cols 73-76: segment identifier
    # cols 77-78: element symbol (right justified)
    
    # Atom name formatting: if 4 chars, left-align; if <4, add space before
    if len(atom_name) == 4:
        atom_field = atom_name
    else:
        atom_field = f" {atom_name:<3}"
    
    line = f"{record:<6}{serial:>5} {atom_field} {res_name:>3} {chain_id:1}{res_num:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}{occupancy:>6.2f}{temp:>6.2f}      {segment:<4}{element:>2}\n"
    
    return line

# Read Phase 2 optimized complex
print("=" * 80)
print("Creating Properly Formatted GLUT1-SDG Complex")
print("=" * 80)

input_file = '/home/pjho3/projects/Drug/scripts/Achive/glut1_tripod_complex.pdb'
output_file = '/home/pjho3/projects/Drug/final_complex/glut1_sdg_complex_optimized_fixed.pdb'

print(f"\nReading: {input_file}")

with open(input_file, 'r') as f:
    lines = f.readlines()

output_lines = []
atom_count = 0
hetatm_count = 0

for line in lines:
    if line.startswith('REMARK') or line.startswith('END'):
        output_lines.append(line)
        continue
    
    if not (line.startswith('ATOM') or line.startswith('HETATM')):
        continue
    
    record = line[0:6].strip()
    
    # Parse fields flexibly
    parts = line.split()
    
    try:
        serial = int(parts[1])
        atom_name = parts[2]
        res_name = parts[3]
        
        # Determine chain ID and residue number
        if parts[4].isalpha() and len(parts[4]) == 1:
            chain_id = parts[4]
            res_num = int(parts[5])
            coord_idx = 6
        else:
            # No chain ID in original - assign 'L' for ligand (SDG)
            if res_name == 'SDG':
                chain_id = 'L'
            else:
                chain_id = ' '
            res_num = int(parts[4])
            coord_idx = 5
        
        x = float(parts[coord_idx])
        y = float(parts[coord_idx + 1])
        z = float(parts[coord_idx + 2])
        
        occupancy = float(parts[coord_idx + 3]) if len(parts) > coord_idx + 3 else 1.00
        temp = float(parts[coord_idx + 4]) if len(parts) > coord_idx + 4 else 0.00
        
        segment = parts[coord_idx + 5] if len(parts) > coord_idx + 5 else ''
        element = parts[coord_idx + 6] if len(parts) > coord_idx + 6 else atom_name[0]
        
        # Format and write
        formatted_line = format_pdb_line(record, serial, atom_name, res_name, chain_id, 
                                        res_num, x, y, z, occupancy, temp, segment, element)
        output_lines.append(formatted_line)
        
        if record == 'ATOM':
            atom_count += 1
        else:
            hetatm_count += 1
            
    except Exception as e:
        print(f"Warning: Skipping line: {line.strip()}")
        print(f"  Error: {e}")
        continue

# Add END record if not present
if not output_lines[-1].startswith('END'):
    output_lines.append('END\n')

print(f"\nWriting: {output_file}")
print(f"  ATOM records: {atom_count}")
print(f"  HETATM records: {hetatm_count}")

with open(output_file, 'w') as f:
    f.writelines(output_lines)

print("\n" + "=" * 80)
print("Verification")
print("=" * 80)

# Show samples
atom_samples = [l for l in output_lines if l.startswith('ATOM')][:3]
hetatm_samples = [l for l in output_lines if l.startswith('HETATM')][:3]

print("\nSample ATOM lines:")
for line in atom_samples:
    print(line.rstrip())

print("\nSample HETATM lines:")
for line in hetatm_samples:
    print(line.rstrip())

print("\n" + "=" * 80)
print("✅ Complex created successfully!")
print("=" * 80)
print(f"\nOutput: {output_file}")
