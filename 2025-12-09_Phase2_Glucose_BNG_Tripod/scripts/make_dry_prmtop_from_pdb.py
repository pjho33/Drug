#!/usr/bin/env python3
"""
Create Amber prmtop from dry_complex.pdb (7432 atoms)
This matches the dry_complex.dcd created by MDAnalysis
"""

import subprocess
import os
import parmed as pmd

print("="*70)
print("Creating Amber prmtop from dry PDB (7432 atoms)")
print("="*70)

# Input
pdb_file = "dry_complex.pdb"
output_prefix = "dry_complex"

print(f"\n📂 Input: {pdb_file}")

# Check atom count
print(f"\n🔍 Checking PDB...")
import MDAnalysis as mda
u = mda.Universe(pdb_file)
print(f"✅ PDB atoms: {len(u.atoms):,}")

# Create tleap script
# Use GAFF2 for SDG ligand, ff14SB for protein, GLYCAM for glycans
tleap_script = f"""
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.GLYCAM_06j-1

# Load structure
mol = loadpdb {pdb_file}

# Check
check mol
charge mol

# Save Amber topology
saveamberparm mol {output_prefix}.parm7 {output_prefix}.rst7

quit
"""

script_file = "tleap_dry.in"
with open(script_file, "w") as f:
    f.write(tleap_script)

print(f"\n📝 Created tleap script: {script_file}")
print(f"\n🔧 Running tleap...")

result = subprocess.run(
    ["tleap", "-f", script_file],
    capture_output=True,
    text=True
)

print(result.stdout)

if result.returncode != 0 or "ERROR" in result.stdout.upper():
    print(f"\n⚠️  tleap warnings/errors detected")
    print(result.stderr)

# Verify output
if os.path.exists(f"{output_prefix}.parm7"):
    parm = pmd.load_file(f"{output_prefix}.parm7")
    print(f"\n✅ Created prmtop: {len(parm.atoms):,} atoms")
    
    if len(parm.atoms) != len(u.atoms):
        print(f"⚠️  Atom count mismatch!")
        print(f"   PDB: {len(u.atoms)}")
        print(f"   Prmtop: {len(parm.atoms)}")
else:
    print(f"\n❌ Failed to create {output_prefix}.parm7")
    exit(1)

print(f"\n{'='*70}")
print("Complete!")
print(f"{'='*70}")
