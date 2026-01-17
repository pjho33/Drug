#!/usr/bin/env python3
"""
Step 1: Create dry trajectory (protein + glycan + SDG only)
Remove water, ions, lipids for MM-GBSA calculation
"""

import MDAnalysis as mda
import os

print("="*70)
print("Step 1: Creating Dry Trajectory")
print("="*70)

# Input files
TRAJ_DIR = "/home/pjho3tr/Downloads/GlycatedGLUT1SDGComplex260107/openmm"
top = os.path.join(TRAJ_DIR, "step5_input.pdb")
traj = os.path.join(TRAJ_DIR, "production.dcd")

# Output directory
OUTPUT_DIR = "/home/pjho3tr/projects/Drug/mmpbsa"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print(f"\n📂 Loading trajectory...")
u = mda.Universe(top, traj)
print(f"✅ Loaded: {len(u.atoms):,} atoms, {len(u.trajectory):,} frames")

# Glycan residue names (common glycosylation residues)
glycan_resnames = [
    "NAG", "BGLCNA",  # N-acetylglucosamine
    "BMA", "MAN", "AMAN", "BMAN",  # Mannose
    "GAL", "BGAL",  # Galactose
    "FUC", "AFUC",  # Fucose
    "SIA", "NEU5AC", "ANEUA",  # Sialic acid
    "NDG", "BGLC", "AGLC",  # Glucose variants
    "A2G", "A2M"  # Other common glycans
]

# Build selection string
print(f"\n🔍 Building selection...")
sel_parts = ["protein", "resname SDG"]
for resname in glycan_resnames:
    sel_parts.append(f"resname {resname}")

sel_string = " or ".join(sel_parts)
print(f"  Selection: {sel_string[:100]}...")

# Select atoms
ag = u.select_atoms(sel_string)
print(f"\n✅ Selected atoms: {len(ag):,}")
print(f"  Protein atoms: {len(u.select_atoms('protein')):,}")
print(f"  Ligand (SDG) atoms: {len(u.select_atoms('resname SDG')):,}")

# Check for glycan
glycan_sel = " or ".join([f"resname {r}" for r in glycan_resnames])
glycan_atoms = u.select_atoms(glycan_sel)
print(f"  Glycan atoms: {len(glycan_atoms):,}")

if len(ag) == 0:
    print("\n❌ Error: No atoms selected!")
    exit(1)

# Write dry complex PDB
output_pdb = os.path.join(OUTPUT_DIR, "dry_complex.pdb")
print(f"\n📝 Writing dry complex PDB...")
ag.write(output_pdb)
print(f"✅ Wrote: {output_pdb}")

# Write dry complex DCD (all frames)
output_dcd = os.path.join(OUTPUT_DIR, "dry_complex.dcd")
print(f"\n📝 Writing dry complex DCD (10,000 frames)...")
print("  This may take several minutes...")

with mda.Writer(output_dcd, ag.n_atoms) as w:
    for i, ts in enumerate(u.trajectory):
        if i % 1000 == 0:
            print(f"    Frame {i:,}/{len(u.trajectory):,} ({100*i/len(u.trajectory):.1f}%)")
        w.write(ag)

print(f"✅ Wrote: {output_dcd}")

# Verify output
u_dry = mda.Universe(output_pdb, output_dcd)
print(f"\n✅ Verification:")
print(f"  Atoms: {len(u_dry.atoms):,}")
print(f"  Frames: {len(u_dry.trajectory):,}")

print(f"\n{'='*70}")
print("Step 1 Complete!")
print(f"{'='*70}")
print(f"\n📁 Output files:")
print(f"  {output_pdb}")
print(f"  {output_dcd}")
