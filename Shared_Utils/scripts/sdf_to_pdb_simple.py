#!/usr/bin/env python3
"""
Simple SDF to PDB converter without external dependencies
"""

def parse_sdf(filename):
    """Parse SDF file and extract atom coordinates"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Find counts line (line 4 in SDF format)
    counts_line = lines[3].strip()
    num_atoms = int(counts_line.split()[0])
    num_bonds = int(counts_line.split()[1])
    
    # Parse atom block (starts at line 5)
    atoms = []
    for i in range(4, 4 + num_atoms):
        parts = lines[i].split()
        x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
        element = parts[3]
        atoms.append((x, y, z, element))
    
    return atoms

def write_pdb(atoms, filename):
    """Write atoms to PDB format"""
    with open(filename, 'w') as f:
        f.write("REMARK   Tripod molecule - (L-glucose-PEG6)3-TRIS\n")
        f.write("REMARK   Generated from Tripod.sdf\n")
        
        for i, (x, y, z, element) in enumerate(atoms, 1):
            # PDB format: ATOM serial name resName chainID resSeq x y z occupancy tempFactor element
            f.write(f"HETATM{i:5d}  {element:<3s} TRP A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n")
        
        f.write("END\n")

# Main conversion
try:
    atoms = parse_sdf('Tripod.sdf')
    write_pdb(atoms, 'tripod_only.pdb')
    print(f"Successfully converted Tripod.sdf to tripod_only.pdb")
    print(f"Number of atoms: {len(atoms)}")
    
    # Count elements
    elements = {}
    for atom in atoms:
        elem = atom[3]
        elements[elem] = elements.get(elem, 0) + 1
    
    print("Element composition:")
    for elem, count in sorted(elements.items()):
        print(f"  {elem}: {count}")
        
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
