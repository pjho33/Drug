#!/usr/bin/env python3
"""
Convert MDAnalysis DCD to Amber NetCDF format
This ensures topology and trajectory match perfectly
"""

import MDAnalysis as mda
import os

print("="*70)
print("Converting DCD to Amber NetCDF format")
print("="*70)

# Input files (MDAnalysis created)
pdb_file = "dry_complex.pdb"
dcd_file = "dry_complex.dcd"

# Output
nc_file = "dry_complex.nc"

print(f"\n📂 Loading trajectory...")
u = mda.Universe(pdb_file, dcd_file)
print(f"✅ Loaded: {len(u.atoms):,} atoms, {len(u.trajectory):,} frames")

print(f"\n📝 Writing NetCDF trajectory...")
print("  This may take several minutes...")

with mda.Writer(nc_file, u.atoms.n_atoms) as w:
    for i, ts in enumerate(u.trajectory):
        if i % 1000 == 0:
            print(f"    Frame {i:,}/{len(u.trajectory):,} ({100*i/len(u.trajectory):.1f}%)")
        w.write(u.atoms)

print(f"✅ Wrote: {nc_file}")

# Verify
u_nc = mda.Universe(pdb_file, nc_file)
print(f"\n✅ Verification:")
print(f"  Atoms: {len(u_nc.atoms):,}")
print(f"  Frames: {len(u_nc.trajectory):,}")

print(f"\n{'='*70}")
print("Complete!")
print(f"{'='*70}")
