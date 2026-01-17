#!/usr/bin/env python3
"""
Fix PDB format - ensure chain ID at column 22
"""

import sys

def fix_pdb_file(input_file, output_file):
    """Fix PDB file format"""
    
    print(f"Processing: {input_file}")
    print(f"Output: {output_file}")
    
    atom_count = 0
    hetatm_count = 0
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Pass through non-ATOM/HETATM lines
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                outfile.write(line)
                continue
            
            # Parse by exact column positions (PDB standard)
            try:
                record = line[0:6].strip()
                serial = int(line[6:11].strip())
                name = line[12:16].strip()
                altLoc = line[16:17].strip() if len(line) > 16 else ''
                resName = line[17:20].strip()
                chainID = line[21:22].strip() if len(line) > 21 else ''
                resSeq = int(line[22:26].strip())
                iCode = line[26:27].strip() if len(line) > 26 else ''
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                occupancy = float(line[54:60].strip()) if len(line) > 54 and line[54:60].strip() else 1.00
                tempFactor = float(line[60:66].strip()) if len(line) > 60 and line[60:66].strip() else 0.00
                element = line[76:78].strip() if len(line) > 76 else ''
                
                # Default chain ID if missing
                if not chainID:
                    if resName in ['TRP', 'SDG', 'UNL', 'LIG']:
                        chainID = 'L'
                    else:
                        chainID = 'A'
                
                # Write in standard PDB format
                # Key: Chain ID MUST be at column 22 (position 21 in 0-indexed)
                new_line = (
                    f"{record:<6s}"      # 1-6
                    f"{serial:>5d}"      # 7-11
                    f" "                 # 12
                    f"{name:<4s}"        # 13-16
                    f"{altLoc:1s}"       # 17
                    f"{resName:>3s}"     # 18-20
                    f" "                 # 21
                    f"{chainID:1s}"      # 22 ← THIS IS THE KEY COLUMN!
                    f"{resSeq:>4d}"      # 23-26
                    f"{iCode:1s}"        # 27
                    f"   "               # 28-30
                    f"{x:8.3f}"          # 31-38
                    f"{y:8.3f}"          # 39-46
                    f"{z:8.3f}"          # 47-54
                    f"{occupancy:6.2f}"  # 55-60
                    f"{tempFactor:6.2f}" # 61-66
                    f"          "        # 67-76
                    f"{element:>2s}"     # 77-78
                    f"\n"
                )
                
                outfile.write(new_line)
                
                if record == 'ATOM':
                    atom_count += 1
                elif record == 'HETATM':
                    hetatm_count += 1
                    
            except Exception as e:
                print(f"Warning: Could not parse line: {line.strip()}")
                print(f"  Error: {e}")
                outfile.write(line)  # Write original line if parsing fails
    
    print(f"  ATOM records: {atom_count}")
    print(f"  HETATM records: {hetatm_count}")
    print("✅ Done!\n")

def main():
    print("=" * 80)
    print("PDB Format Fixer - Chain ID at Column 22")
    print("=" * 80)
    print()
    
    files = [
        ("/home/pjho3/projects/Drug/phase3_charmm_gui_submission/experimental/glut1_tripod_complex.pdb",
         "/home/pjho3/projects/Drug/phase3_charmm_gui_submission/experimental/glut1_tripod_complex_fixed.pdb"),
        ("/home/pjho3/projects/Drug/phase3_charmm_gui_submission/control/glut1_tripod_complex_control.pdb",
         "/home/pjho3/projects/Drug/phase3_charmm_gui_submission/control/glut1_tripod_complex_control_fixed.pdb"),
    ]
    
    import os
    for input_file, output_file in files:
        if os.path.exists(input_file):
            fix_pdb_file(input_file, output_file)
        else:
            print(f"⚠️  File not found: {input_file}\n")
    
    print("=" * 80)
    print("✅ All files processed!")
    print("=" * 80)
    print("\nFixed files ready for CHARMM-GUI:")
    for _, output_file in files:
        if os.path.exists(output_file):
            print(f"  ✓ {output_file}")
    
    print("\n" + "=" * 80)

if __name__ == "__main__":
    main()
