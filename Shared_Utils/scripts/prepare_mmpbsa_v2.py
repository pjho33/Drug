#!/usr/bin/env python3
"""
Prepare files for MMPBSA.py calculation with proper ligand parameterization
"""

import MDAnalysis as mda
from pathlib import Path
import subprocess
import sys

# Paths
BASE_DIR = Path("/home/pjho3/projects/Drug/final_complex/controlcomplex/openmm")
MMPBSA_DIR = Path("/home/pjho3/projects/Drug/final_complex/mmpbsa")
MMPBSA_DIR.mkdir(exist_ok=True)

print("=" * 80)
print("Preparing Files for MMPBSA.py (with ligand parameterization)")
print("=" * 80)
print()

# Input files
PSF_FILE = BASE_DIR / "step5_input.psf"
DCD_FILE = BASE_DIR / "step7_1.dcd"

print(f"Input files:")
print(f"  PSF: {PSF_FILE}")
print(f"  DCD: {DCD_FILE}")
print()

# Load universe
print("Loading trajectory...")
u = mda.Universe(str(PSF_FILE), str(DCD_FILE))
print(f"✅ Loaded: {len(u.trajectory)} frames")
print()

# Identify components
protein = u.select_atoms("protein")
ligand = u.select_atoms("resname SDG")
if len(ligand) == 0:
    ligand = u.select_atoms("not protein and not resname TIP3 and not resname SOD and not resname CLA and not resname POT")

print(f"System components:")
print(f"  Protein: {len(protein)} atoms")
print(f"  Ligand: {len(ligand)} atoms ({ligand.resnames[0]})")
print()

# Get first frame
u.trajectory[0]

# ============================================================================
# Step 1: Create PDB files
# ============================================================================
print("=" * 80)
print("Step 1: Creating PDB files")
print("=" * 80)
print()

complex_sel = u.select_atoms("protein or resname SDG")

print("Writing PDB files...")
complex_sel.write(str(MMPBSA_DIR / "complex.pdb"))
protein.write(str(MMPBSA_DIR / "receptor.pdb"))
ligand.write(str(MMPBSA_DIR / "ligand.pdb"))
print(f"  ✅ complex.pdb")
print(f"  ✅ receptor.pdb")
print(f"  ✅ ligand.pdb")
print()

# ============================================================================
# Step 2: Parameterize ligand with antechamber
# ============================================================================
print("=" * 80)
print("Step 2: Parameterizing ligand with antechamber")
print("=" * 80)
print()

print("Running antechamber to generate GAFF parameters...")
antechamber_cmd = [
    "antechamber",
    "-i", str(MMPBSA_DIR / "ligand.pdb"),
    "-fi", "pdb",
    "-o", str(MMPBSA_DIR / "ligand.mol2"),
    "-fo", "mol2",
    "-c", "bcc",  # AM1-BCC charges
    "-s", "2",    # Verbosity
    "-nc", "0",   # Net charge (adjust if needed)
    "-at", "gaff2"
]

result = subprocess.run(antechamber_cmd, capture_output=True, text=True, cwd=str(MMPBSA_DIR))
if result.returncode == 0:
    print(f"  ✅ ligand.mol2 created")
else:
    print(f"  ⚠️  antechamber warning:")
    print(result.stdout[-500:] if len(result.stdout) > 500 else result.stdout)
print()

print("Running parmchk2 to check parameters...")
parmchk_cmd = [
    "parmchk2",
    "-i", str(MMPBSA_DIR / "ligand.mol2"),
    "-f", "mol2",
    "-o", str(MMPBSA_DIR / "ligand.frcmod"),
    "-s", "gaff2"
]

result = subprocess.run(parmchk_cmd, capture_output=True, text=True, cwd=str(MMPBSA_DIR))
if result.returncode == 0:
    print(f"  ✅ ligand.frcmod created")
else:
    print(f"  ⚠️  parmchk2 warning:")
    print(result.stdout)
print()

# ============================================================================
# Step 3: Create topology files with tleap
# ============================================================================
print("=" * 80)
print("Step 3: Creating Amber topology files")
print("=" * 80)
print()

# Complex topology
tleap_complex = f"""source leaprc.protein.ff14SB
source leaprc.gaff2

# Load ligand parameters
loadamberparams {MMPBSA_DIR / "ligand.frcmod"}
LIG = loadmol2 {MMPBSA_DIR / "ligand.mol2"}

# Load receptor
REC = loadpdb {MMPBSA_DIR / "receptor.pdb"}

# Create complex
complex = combine {{REC LIG}}

# Save topology
saveamberparm complex {MMPBSA_DIR / "complex.prmtop"} {MMPBSA_DIR / "complex.inpcrd"}
quit
"""

with open(MMPBSA_DIR / "tleap_complex.in", 'w') as f:
    f.write(tleap_complex)

print("Creating complex topology...")
result = subprocess.run(
    ["tleap", "-f", "tleap_complex.in"],
    cwd=str(MMPBSA_DIR),
    capture_output=True,
    text=True
)

if (MMPBSA_DIR / "complex.prmtop").exists():
    print(f"  ✅ complex.prmtop")
else:
    print(f"  ❌ Failed to create complex.prmtop")
    print(result.stdout[-1000:])
print()

# Receptor topology
tleap_receptor = f"""source leaprc.protein.ff14SB

REC = loadpdb {MMPBSA_DIR / "receptor.pdb"}
saveamberparm REC {MMPBSA_DIR / "receptor.prmtop"} {MMPBSA_DIR / "receptor.inpcrd"}
quit
"""

with open(MMPBSA_DIR / "tleap_receptor.in", 'w') as f:
    f.write(tleap_receptor)

print("Creating receptor topology...")
result = subprocess.run(
    ["tleap", "-f", "tleap_receptor.in"],
    cwd=str(MMPBSA_DIR),
    capture_output=True,
    text=True
)

if (MMPBSA_DIR / "receptor.prmtop").exists():
    print(f"  ✅ receptor.prmtop")
else:
    print(f"  ❌ Failed to create receptor.prmtop")
print()

# Ligand topology
tleap_ligand = f"""source leaprc.gaff2

loadamberparams {MMPBSA_DIR / "ligand.frcmod"}
LIG = loadmol2 {MMPBSA_DIR / "ligand.mol2"}
saveamberparm LIG {MMPBSA_DIR / "ligand.prmtop"} {MMPBSA_DIR / "ligand.inpcrd"}
quit
"""

with open(MMPBSA_DIR / "tleap_ligand.in", 'w') as f:
    f.write(tleap_ligand)

print("Creating ligand topology...")
result = subprocess.run(
    ["tleap", "-f", "tleap_ligand.in"],
    cwd=str(MMPBSA_DIR),
    capture_output=True,
    text=True
)

if (MMPBSA_DIR / "ligand.prmtop").exists():
    print(f"  ✅ ligand.prmtop")
else:
    print(f"  ❌ Failed to create ligand.prmtop")
print()

# ============================================================================
# Step 4: Convert trajectory
# ============================================================================
print("=" * 80)
print("Step 4: Converting trajectory to Amber format")
print("=" * 80)
print()

print("Writing trajectory (complex, no solvent)...")
with mda.Writer(str(MMPBSA_DIR / "trajectory.nc"), complex_sel.n_atoms) as W:
    for ts in u.trajectory:
        W.write(complex_sel)
print(f"  ✅ trajectory.nc ({len(u.trajectory)} frames)")
print()

# ============================================================================
# Step 5: Create MMPBSA input
# ============================================================================
print("=" * 80)
print("Step 5: Creating MMPBSA.py input file")
print("=" * 80)
print()

mmpbsa_input = """&general
startframe=1
endframe=10
interval=1
verbose=2
/

&gb
igb=5, saltcon=0.150
/

&pb
istrng=0.150, fillratio=4.0
/
"""

with open(MMPBSA_DIR / "mmpbsa.in", 'w') as f:
    f.write(mmpbsa_input)

print(f"✅ mmpbsa.in created")
print()

# ============================================================================
# Summary
# ============================================================================
print("=" * 80)
print("✅ Preparation Complete!")
print("=" * 80)
print()
print(f"Output directory: {MMPBSA_DIR}")
print()
print("Files created:")
print("  ✅ complex.pdb, receptor.pdb, ligand.pdb")
print("  ✅ ligand.mol2, ligand.frcmod (GAFF2 parameters)")
print("  ✅ complex.prmtop, receptor.prmtop, ligand.prmtop")
print("  ✅ trajectory.nc")
print("  ✅ mmpbsa.in")
print()
print("To run MMPBSA.py:")
print(f"  cd {MMPBSA_DIR}")
print("  MMPBSA.py -O -i mmpbsa.in -o mmpbsa_results.dat \\")
print("    -sp complex.prmtop -cp complex.prmtop \\")
print("    -rp receptor.prmtop -lp ligand.prmtop \\")
print("    -y trajectory.nc")
print()
