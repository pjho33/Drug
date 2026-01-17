#!/usr/bin/env python3
"""
Prepare files for MMPBSA.py calculation
Convert CHARMM trajectory to Amber format and create topology files
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
print("Preparing Files for MMPBSA.py")
print("=" * 80)
print()

# Input files
PSF_FILE = BASE_DIR / "step5_input.psf"
PDB_FILE = BASE_DIR / "step5_input.pdb"
DCD_FILE = BASE_DIR / "step7_1.dcd"

print(f"Input files:")
print(f"  PSF: {PSF_FILE}")
print(f"  PDB: {PDB_FILE}")
print(f"  DCD: {DCD_FILE}")
print()

# Load universe
print("Loading trajectory...")
u = mda.Universe(str(PSF_FILE), str(DCD_FILE))
print(f"✅ Loaded: {len(u.trajectory)} frames")
print(f"   Atoms: {len(u.atoms)}")
print()

# Identify components
print("Identifying system components...")
protein = u.select_atoms("protein")
print(f"  Protein atoms: {len(protein)}")

# Find ligand
ligand = u.select_atoms("resname SDG")
if len(ligand) == 0:
    ligand = u.select_atoms("not protein and not resname TIP3 and not resname SOD and not resname CLA and not resname POT")
print(f"  Ligand atoms: {len(ligand)}")
print(f"  Ligand resname: {ligand.resnames[0] if len(ligand) > 0 else 'Not found'}")

# Complex (protein + ligand, no water/ions)
complex_sel = u.select_atoms("protein or resname SDG or (not protein and not resname TIP3 and not resname SOD and not resname CLA and not resname POT)")
print(f"  Complex atoms (no solvent): {len(complex_sel)}")
print()

# ============================================================================
# Create PDB files for complex, receptor, and ligand
# ============================================================================
print("=" * 80)
print("Creating PDB files for MMPBSA")
print("=" * 80)
print()

# Get first frame
u.trajectory[0]

# Complex PDB (protein + ligand)
print("Writing complex.pdb...")
complex_sel.write(str(MMPBSA_DIR / "complex.pdb"))
print(f"  ✅ {MMPBSA_DIR / 'complex.pdb'}")

# Receptor PDB (protein only)
print("Writing receptor.pdb...")
protein.write(str(MMPBSA_DIR / "receptor.pdb"))
print(f"  ✅ {MMPBSA_DIR / 'receptor.pdb'}")

# Ligand PDB (ligand only)
print("Writing ligand.pdb...")
ligand.write(str(MMPBSA_DIR / "ligand.pdb"))
print(f"  ✅ {MMPBSA_DIR / 'ligand.pdb'}")
print()

# ============================================================================
# Convert trajectory to Amber format (NetCDF)
# ============================================================================
print("=" * 80)
print("Converting trajectory to Amber NetCDF format")
print("=" * 80)
print()

print("Writing trajectory (complex, no solvent)...")
with mda.Writer(str(MMPBSA_DIR / "trajectory.nc"), complex_sel.n_atoms) as W:
    for ts in u.trajectory:
        W.write(complex_sel)
print(f"  ✅ {MMPBSA_DIR / 'trajectory.nc'}")
print(f"  Frames: {len(u.trajectory)}")
print()

# ============================================================================
# Create topology files using tleap
# ============================================================================
print("=" * 80)
print("Creating Amber topology files with tleap")
print("=" * 80)
print()

# Create tleap input for complex
tleap_complex = f"""source leaprc.protein.ff14SB
source leaprc.gaff2

# Load PDBs
complex = loadpdb {MMPBSA_DIR / "complex.pdb"}

# Save topology and coordinates
saveamberparm complex {MMPBSA_DIR / "complex.prmtop"} {MMPBSA_DIR / "complex.inpcrd"}
quit
"""

tleap_input_file = MMPBSA_DIR / "tleap_complex.in"
with open(tleap_input_file, 'w') as f:
    f.write(tleap_complex)

print("Running tleap for complex...")
print(f"  Input: {tleap_input_file}")
result = subprocess.run(
    ["tleap", "-f", str(tleap_input_file)],
    cwd=str(MMPBSA_DIR),
    capture_output=True,
    text=True
)

if result.returncode == 0:
    print(f"  ✅ {MMPBSA_DIR / 'complex.prmtop'}")
    print(f"  ✅ {MMPBSA_DIR / 'complex.inpcrd'}")
else:
    print(f"  ⚠️  tleap warning/error:")
    print(result.stdout)
    print(result.stderr)
print()

# Create tleap input for receptor
tleap_receptor = f"""source leaprc.protein.ff14SB

# Load PDB
receptor = loadpdb {MMPBSA_DIR / "receptor.pdb"}

# Save topology and coordinates
saveamberparm receptor {MMPBSA_DIR / "receptor.prmtop"} {MMPBSA_DIR / "receptor.inpcrd"}
quit
"""

tleap_input_file = MMPBSA_DIR / "tleap_receptor.in"
with open(tleap_input_file, 'w') as f:
    f.write(tleap_receptor)

print("Running tleap for receptor...")
print(f"  Input: {tleap_input_file}")
result = subprocess.run(
    ["tleap", "-f", str(tleap_input_file)],
    cwd=str(MMPBSA_DIR),
    capture_output=True,
    text=True
)

if result.returncode == 0:
    print(f"  ✅ {MMPBSA_DIR / 'receptor.prmtop'}")
    print(f"  ✅ {MMPBSA_DIR / 'receptor.inpcrd'}")
else:
    print(f"  ⚠️  tleap warning/error:")
    print(result.stdout)
    print(result.stderr)
print()

# Create tleap input for ligand
tleap_ligand = f"""source leaprc.gaff2

# Load PDB
ligand = loadpdb {MMPBSA_DIR / "ligand.pdb"}

# Save topology and coordinates
saveamberparm ligand {MMPBSA_DIR / "ligand.prmtop"} {MMPBSA_DIR / "ligand.inpcrd"}
quit
"""

tleap_input_file = MMPBSA_DIR / "tleap_ligand.in"
with open(tleap_input_file, 'w') as f:
    f.write(tleap_ligand)

print("Running tleap for ligand...")
print(f"  Input: {tleap_input_file}")
result = subprocess.run(
    ["tleap", "-f", str(tleap_input_file)],
    cwd=str(MMPBSA_DIR),
    capture_output=True,
    text=True
)

if result.returncode == 0:
    print(f"  ✅ {MMPBSA_DIR / 'ligand.prmtop'}")
    print(f"  ✅ {MMPBSA_DIR / 'ligand.inpcrd'}")
else:
    print(f"  ⚠️  tleap warning/error:")
    print(result.stdout)
    print(result.stderr)
print()

# ============================================================================
# Create MMPBSA.py input file
# ============================================================================
print("=" * 80)
print("Creating MMPBSA.py input file")
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

mmpbsa_input_file = MMPBSA_DIR / "mmpbsa.in"
with open(mmpbsa_input_file, 'w') as f:
    f.write(mmpbsa_input)

print(f"✅ MMPBSA input file created: {mmpbsa_input_file}")
print()
print("Input file contents:")
print(mmpbsa_input)

# ============================================================================
# Summary
# ============================================================================
print("=" * 80)
print("Preparation Complete!")
print("=" * 80)
print()
print(f"All files saved to: {MMPBSA_DIR}")
print()
print("Files created:")
print("  PDB files:")
print("    - complex.pdb")
print("    - receptor.pdb")
print("    - ligand.pdb")
print("  Topology files:")
print("    - complex.prmtop")
print("    - receptor.prmtop")
print("    - ligand.prmtop")
print("  Trajectory:")
print("    - trajectory.nc")
print("  MMPBSA input:")
print("    - mmpbsa.in")
print()
print("To run MMPBSA.py:")
print(f"  cd {MMPBSA_DIR}")
print("  MMPBSA.py -O -i mmpbsa.in -o mmpbsa_results.dat \\")
print("    -sp complex.prmtop -cp complex.prmtop \\")
print("    -rp receptor.prmtop -lp ligand.prmtop \\")
print("    -y trajectory.nc")
print()
