#!/usr/bin/env python3
"""
Prepare and run MMPBSA.py with Amber format files from CHARMM-GUI
"""

import MDAnalysis as mda
from pathlib import Path
import subprocess
import sys

# Paths
AMBER_DIR = Path("/home/pjho3/projects/Drug/final_complex/GLUT1SDGComplex260110/amber")
MMPBSA_DIR = Path("/home/pjho3/projects/Drug/final_complex/mmpbsa_amber")
MMPBSA_DIR.mkdir(exist_ok=True)

print("=" * 80)
print("MMPBSA.py with Amber Format Files")
print("=" * 80)
print()

# Input files
PARM7_FILE = AMBER_DIR / "step5_input.parm7"
PDB_FILE = AMBER_DIR / "step5_input.pdb"

print(f"Input files:")
print(f"  Topology: {PARM7_FILE}")
print(f"  PDB: {PDB_FILE}")
print()

# Check if trajectory exists
traj_files = list(AMBER_DIR.glob("*.nc")) + list(AMBER_DIR.glob("*.mdcrd"))
if traj_files:
    print(f"Found trajectory files: {[f.name for f in traj_files]}")
else:
    print("⚠️  No trajectory found. Will use PDB coordinates only.")

print()

# ============================================================================
# Step 1: Load structure and identify components
# ============================================================================
print("=" * 80)
print("Step 1: Identifying system components")
print("=" * 80)
print()

u = mda.Universe(str(PARM7_FILE), str(PDB_FILE))
print(f"Total atoms: {len(u.atoms)}")

protein = u.select_atoms("protein")
print(f"Protein atoms: {len(protein)}")

ligand = u.select_atoms("resname SDG")
print(f"Ligand atoms: {len(ligand)} (resname: {ligand.resnames[0] if len(ligand) > 0 else 'Not found'})")

# Get residue numbers
protein_residues = sorted(set(protein.resids))
ligand_residues = sorted(set(ligand.resids))

print(f"\nProtein residue range: {protein_residues[0]}-{protein_residues[-1]}")
print(f"Ligand residue: {ligand_residues[0] if ligand_residues else 'N/A'}")
print()

# ============================================================================
# Step 2: Create MMPBSA input file
# ============================================================================
print("=" * 80)
print("Step 2: Creating MMPBSA.py input file")
print("=" * 80)
print()

# MMPBSA input with both GB and PB calculations
mmpbsa_input = f"""Input file for running MMPBSA
&general
startframe=1
endframe=1
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
# Step 3: Create receptor and ligand masks
# ============================================================================
print("=" * 80)
print("Step 3: Determining atom masks")
print("=" * 80)
print()

# Amber uses 1-based indexing
receptor_mask = f":{protein_residues[0]}-{protein_residues[-1]}"
ligand_mask = f":{ligand_residues[0]}" if ligand_residues else ":SDG"

print(f"Receptor mask: {receptor_mask}")
print(f"Ligand mask: {ligand_mask}")
print()

# ============================================================================
# Step 4: Run MMPBSA.py
# ============================================================================
print("=" * 80)
print("Step 4: Running MMPBSA.py")
print("=" * 80)
print()

# Use PDB as "trajectory" for single-point calculation
cmd = [
    "MMPBSA.py",
    "-O",  # Overwrite output
    "-i", str(mmpbsa_input_file),
    "-o", str(MMPBSA_DIR / "FINAL_RESULTS_MMPBSA.dat"),
    "-sp", str(PARM7_FILE),  # Solvated topology
    "-cp", str(PARM7_FILE),  # Complex topology
    "-rp", str(PARM7_FILE),  # Receptor topology (will be stripped)
    "-lp", str(PARM7_FILE),  # Ligand topology (will be stripped)
    "-y", str(PDB_FILE),     # Use PDB as trajectory
    "-srp", receptor_mask,   # Receptor mask for stripping
    "-slp", ligand_mask      # Ligand mask for stripping
]

print("Command:")
print(" ".join(cmd))
print()
print("Running MMPBSA.py (this may take several minutes)...")
print("-" * 80)

result = subprocess.run(
    cmd,
    cwd=str(MMPBSA_DIR),
    capture_output=True,
    text=True
)

# Save output
with open(MMPBSA_DIR / "mmpbsa_stdout.log", 'w') as f:
    f.write(result.stdout)
with open(MMPBSA_DIR / "mmpbsa_stderr.log", 'w') as f:
    f.write(result.stderr)

if result.returncode == 0:
    print("✅ MMPBSA.py completed successfully!")
else:
    print(f"❌ MMPBSA.py failed with exit code {result.returncode}")
    print("\nStderr output:")
    print(result.stderr[-1000:])

print()

# ============================================================================
# Step 5: Display results
# ============================================================================
print("=" * 80)
print("Step 5: Results")
print("=" * 80)
print()

results_file = MMPBSA_DIR / "FINAL_RESULTS_MMPBSA.dat"
if results_file.exists():
    print(f"Results saved to: {results_file}")
    print()
    print("=" * 80)
    print("FINAL RESULTS")
    print("=" * 80)
    with open(results_file, 'r') as f:
        print(f.read())
else:
    print("⚠️  Results file not found. Check logs for errors.")
    print(f"\nStdout log: {MMPBSA_DIR / 'mmpbsa_stdout.log'}")
    print(f"Stderr log: {MMPBSA_DIR / 'mmpbsa_stderr.log'}")

print()
print("=" * 80)
print("✅ MMPBSA Analysis Complete!")
print("=" * 80)
print(f"\nAll files saved to: {MMPBSA_DIR}")
print()
