#!/usr/bin/env python3
"""
Add radii to original parm7 and recreate separated topologies
"""

import parmed as pmd
from pathlib import Path
import subprocess

AMBER_DIR = Path("/home/pjho3/projects/Drug/final_complex/GLUT1SDGComplex260110/amber")
MMPBSA_DIR = Path("/home/pjho3/projects/Drug/final_complex/mmpbsa_amber")

print("=" * 80)
print("Adding Radii to Original Topology")
print("=" * 80)
print()

# Load original parm7
original_parm = AMBER_DIR / "step5_input.parm7"
print(f"Loading: {original_parm}")
parm = pmd.load_file(str(original_parm))
print(f"  Atoms: {len(parm.atoms)}")
print(f"  Residues: {len(parm.residues)}")
print()

# Add radii
print("Adding mbondi2 radii...")
parm = pmd.tools.changeRadii(parm, 'mbondi2')
print("  ✅ Radii added")
print()

# Save with radii
parm7_radii = MMPBSA_DIR / "step5_input_radii.parm7"
print(f"Saving to: {parm7_radii}")
parm.save(str(parm7_radii), overwrite=True)
print("  ✅ Saved")
print()

# Now run ante-MMPBSA with radii-enabled topology
print("=" * 80)
print("Running ante-MMPBSA.py with radii-enabled topology")
print("=" * 80)
print()

cmd = [
    "ante-MMPBSA.py",
    "-p", str(parm7_radii),
    "-c", str(MMPBSA_DIR / "complex.prmtop"),
    "-r", str(MMPBSA_DIR / "receptor.prmtop"),
    "-l", str(MMPBSA_DIR / "ligand.prmtop"),
    "-s", ":306",
    "-n", ":1-305"
]

print("Command:")
print(" ".join(cmd))
print()

result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(MMPBSA_DIR))

print(result.stdout)
if result.stderr:
    print("Stderr:", result.stderr)

if result.returncode == 0:
    print()
    print("=" * 80)
    print("✅ Topology files created successfully with radii!")
    print("=" * 80)
    print()
    print("Files created:")
    print(f"  - {MMPBSA_DIR / 'complex.prmtop'}")
    print(f"  - {MMPBSA_DIR / 'receptor.prmtop'}")
    print(f"  - {MMPBSA_DIR / 'ligand.prmtop'}")
else:
    print()
    print("❌ ante-MMPBSA.py failed")
    print(f"Exit code: {result.returncode}")
