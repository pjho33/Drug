#!/usr/bin/env python3
"""
Add implicit radii to topology files for GB calculations
"""

import parmed as pmd
from pathlib import Path

MMPBSA_DIR = Path("/home/pjho3/projects/Drug/final_complex/mmpbsa_amber")

print("=" * 80)
print("Adding Implicit Radii to Topology Files")
print("=" * 80)
print()

# Process each topology file
for topo_name in ['complex', 'receptor', 'ligand']:
    topo_file = MMPBSA_DIR / f"{topo_name}.prmtop"
    
    print(f"Processing {topo_name}.prmtop...")
    
    # Load topology
    parm = pmd.load_file(str(topo_file))
    
    # Add radii (mbondi2 is recommended for GB calculations)
    parm = pmd.tools.addPDB(parm)
    parm = pmd.tools.changeRadii(parm, 'mbondi2')
    
    # Save with radii
    parm.save(str(topo_file), overwrite=True)
    
    print(f"  ✅ Added mbondi2 radii to {topo_name}.prmtop")

print()
print("=" * 80)
print("✅ All topology files updated with implicit radii")
print("=" * 80)
