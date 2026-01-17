#!/usr/bin/env python3
"""
Fix PDB format - ensure all columns are correctly aligned
Especially chain ID at column 22
"""

import sys
import os

def fix_pdb_line(line):
    """Fix a single PDB line to standard format"""
    
    # Skip non-ATOM/HETATM lines
    if not (line.startswith('ATOM') or line.startswith('HETATM')):
        return line
    
    try:
        # Parse the line
        record = line[0:6].strip()  # ATOM or HETATM
        serial = line[6:11].strip()  # Atom serial number
        name = line[12:16].strip()   # Atom name
        altLoc = line[16:17].strip() # Alternate location
        resName = line[17:20].strip() # Residue name
        chainID = line[21:22].strip() # Chain ID
        resSeq = line[22:26].strip()  # Residue sequence number
        iCode = line[26:27].strip()   # Insertion code
        x = line[30:38].strip()       # X coordinate
        y = line[38:46].strip()       # Y coordinate
        z = line[46:54].strip()       # Z coordinate
        occupancy = line[54:60].strip() # Occupancy
        tempFactor = line[60:66].strip() # Temperature factor
        element = line[76:78].strip() if len(line) > 76 else ''  # Element
        
        # Handle missing chain ID - use 'A' as default for ligands
        if not chainID:
            if resName in ['TRP', 'SDG', 'UNL', 'LIG']:
                chainID = 'L'  # Ligand chain
            else:
                chainID = 'A'  # Protein chain
        
        # Format coordinates
        try:
            x_float = float(x)
            y_float = float(y)
            z_float = float(z)
        except:
            return line  # Return original if parsing fails
        
        # Format occupancy and temp factor
        try:
            occ = float(occupancy) if occupancy else 1.00
            temp = float(tempFactor) if tempFactor else 0.00
        except:
            occ = 1.00
            temp = 0.00
        
        # Reconstruct line in standard PDB format
        # Columns:  1-6    7-11  12    13-16 17    18-20 21  22    23-26 27   28-30      31-38      39-46      47-54    55-60  61-66        77-78
        new_line = (
            f"{record:<6s}"           # 1-6: Record name
            f"{int(serial):>5d}"      # 7-11: Atom serial
            f" "                      # 12: Space
            f"{name:<4s}"             # 13-16: Atom name (left-aligned for 1-2 char, right for 4)
            f"{altLoc:1s}"            # 17: Alternate location
            f"{resName:>3s}"          # 18-20: Residue name (right-aligned)
            f" "                      # 21: Space
            f"{chainID:1s}"           # 22: Chain ID (THIS IS THE KEY!)
            f"{int(resSeq):>4d}"      # 23-26: Residue sequence number
            f"{iCode:1s}"             # 27: Insertion code
            f"   "                    # 28-30: Spaces
            f"{x_float:8.3f}"         # 31-38: X coordinate
            f"{y_float:8.3f}"         # 39-46: Y coordinate
            f"{z_float:8.3f}"         # 47-54: Z coordinate
            f"{occ:6.2f}"             # 55-60: Occupancy
            f"{temp:6.2f}"            # 61-66: Temperature factor
            f"          "             # 67-76: Spaces
            f"{element:>2s}"          # 77-78: Element symbol
            f"\n"
        )
        
        return new_line
        
    except Exception as e:
        print(f"Warning: Could not parse line: {line.strip()}")
        print(f"  Error: {e}")
        return line

def fix_pdb_file(input_file, output_file):
    """Fix entire PDB file"""
    
    print(f"Processing: {input_file}")
    print(f"Output to: {output_file}")
    
    atom_count = 0
    hetatm_count = 0
    other_count = 0
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('ATOM'):
                atom_count += 1
                fixed_line = fix_pdb_line(line)
                outfile.write(fixed_line)
            elif line.startswith('HETATM'):
                hetatm_count += 1
                fixed_line = fix_pdb_line(line)
                outfile.write(fixed_line)
            else:
                other_count += 1
                outfile.write(line)
    
    print(f"  ATOM records: {atom_count}")
    print(f"  HETATM records: {hetatm_count}")
    print(f"  Other lines: {other_count}")
    print(f"✅ Done!\n")

def main():
    print("=" * 80)
    print("PDB Format Fixer - Column Alignment")
    print("=" * 80)
    print("\nThis script fixes PDB files to ensure:")
    print("  - Chain ID is at column 22 (REQUIRED by CHARMM-GUI)")
    print("  - All coordinates are properly aligned")
    print("  - Standard PDB format compliance")
    print()
    
    # File paths
    files_to_fix = [
        ("/home/pjho3/projects/Drug/phase3_charmm_gui_submission/experimental/glut1_tripod_complex.pdb",
         "/home/pjho3/projects/Drug/phase3_charmm_gui_submission/experimental/glut1_tripod_complex_fixed.pdb"),
        ("/home/pjho3/projects/Drug/phase3_charmm_gui_submission/control/glut1_tripod_complex_control.pdb",
         "/home/pjho3/projects/Drug/phase3_charmm_gui_submission/control/glut1_tripod_complex_control_fixed.pdb"),
    ]
    
    for input_file, output_file in files_to_fix:
        if os.path.exists(input_file):
            fix_pdb_file(input_file, output_file)
        else:
            print(f"⚠️  File not found: {input_file}")
            print()
    
    print("=" * 80)
    print("✅ All files processed!")
    print("=" * 80)
    print("\nFixed files:")
    for _, output_file in files_to_fix:
        if os.path.exists(output_file):
            print(f"  ✓ {output_file}")
    
    print("\n📋 Next step:")
    print("  Upload the *_fixed.pdb files to CHARMM-GUI")
    print("=" * 80)

if __name__ == "__main__":
    main()
