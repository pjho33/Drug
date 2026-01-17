#!/usr/bin/env python3
"""
Additional Trajectory Analysis: RMSF and Hydrogen Bonds
"""

import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSF
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
import numpy as np
import matplotlib.pyplot as plt
import os

print("="*70)
print("Additional Analysis: RMSF & Hydrogen Bonds")
print("="*70)

# Paths
TRAJ_DIR = "/home/pjho3tr/Downloads/GlycatedGLUT1SDGComplex260107/openmm"
PSF_FILE = os.path.join(TRAJ_DIR, "step5_input.psf")
DCD_FILE = os.path.join(TRAJ_DIR, "production.dcd")
OUTPUT_DIR = "/home/pjho3tr/projects/Drug/analysis"

print(f"\n📂 Loading trajectory...")
u = mda.Universe(PSF_FILE, DCD_FILE)
print(f"✅ Loaded: {len(u.trajectory):,} frames")

protein = u.select_atoms("protein")
protein_ca = u.select_atoms("protein and name CA")
ligand = u.select_atoms("resname SDG")

# ============================================================================
# 1. RMSF Analysis (Per-Residue Flexibility)
# ============================================================================
print(f"\n{'='*70}")
print("1. RMSF Analysis")
print(f"{'='*70}")

print("  Calculating per-residue RMSF...")
rmsf_analysis = RMSF(protein_ca).run()

# Save RMSF data
rmsf_data = np.column_stack([
    protein_ca.resids,
    rmsf_analysis.results.rmsf
])
np.savetxt(
    os.path.join(OUTPUT_DIR, "rmsf_data.txt"),
    rmsf_data,
    header="Residue RMSF(Å)",
    fmt="%d %.4f"
)

# Plot RMSF
plt.figure(figsize=(14, 6))
plt.plot(rmsf_data[:, 0], rmsf_data[:, 1], 'g-', linewidth=1.5)
plt.xlabel('Residue Number', fontsize=12)
plt.ylabel('RMSF (Å)', fontsize=12)
plt.title('Per-Residue Flexibility (RMSF)', fontsize=14, fontweight='bold')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "rmsf_plot.png"), dpi=300)
plt.close()
print(f"✅ RMSF plot saved")

print(f"\n📊 RMSF Statistics:")
print(f"  Mean RMSF: {np.mean(rmsf_data[:, 1]):.2f} Å")
print(f"  Max RMSF: {np.max(rmsf_data[:, 1]):.2f} Å (Residue {rmsf_data[np.argmax(rmsf_data[:, 1]), 0]:.0f})")

# Identify highly flexible residues (RMSF > 2.0 Å)
flexible_residues = rmsf_data[rmsf_data[:, 1] > 2.0]
if len(flexible_residues) > 0:
    print(f"\n  Highly flexible residues (RMSF > 2.0 Å):")
    for res_id, res_rmsf in flexible_residues[:10]:  # Show top 10
        print(f"    Residue {res_id:.0f}: {res_rmsf:.2f} Å")

# ============================================================================
# 2. Hydrogen Bond Analysis
# ============================================================================
print(f"\n{'='*70}")
print("2. Hydrogen Bond Analysis")
print(f"{'='*70}")

print("  Analyzing protein-ligand hydrogen bonds...")
print("  (This may take several minutes for 10,000 frames)")

hbond_analysis = HydrogenBondAnalysis(
    universe=u,
    donors_sel="protein",
    hydrogens_sel="protein",
    acceptors_sel="resname SDG",
    d_a_cutoff=3.5,
    d_h_a_angle_cutoff=150
)
hbond_analysis.run(verbose=True)

# Count hydrogen bonds per frame
hbond_counts = hbond_analysis.count_by_time()
time_ns = np.arange(len(hbond_counts)) * u.trajectory.dt / 1000

# Save hydrogen bond data
hbond_data = np.column_stack([time_ns, hbond_counts])
np.savetxt(
    os.path.join(OUTPUT_DIR, "hbond_data.txt"),
    hbond_data,
    header="Time(ns) HBond_Count",
    fmt="%.4f %d"
)

# Plot hydrogen bonds
plt.figure(figsize=(12, 6))
plt.plot(time_ns, hbond_counts, 'purple', linewidth=1.5)
plt.xlabel('Time (ns)', fontsize=12)
plt.ylabel('Number of Hydrogen Bonds', fontsize=12)
plt.title('Protein-Ligand Hydrogen Bonds', fontsize=14, fontweight='bold')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "hbond_plot.png"), dpi=300)
plt.close()
print(f"✅ Hydrogen bond plot saved")

print(f"\n📊 Hydrogen Bond Statistics:")
print(f"  Average H-bonds: {np.mean(hbond_counts):.2f} ± {np.std(hbond_counts):.2f}")
print(f"  Max H-bonds: {np.max(hbond_counts)}")
print(f"  Min H-bonds: {np.min(hbond_counts)}")
print(f"  Occupancy (>0 H-bonds): {100 * np.sum(hbond_counts > 0) / len(hbond_counts):.1f}%")

# ============================================================================
# Summary
# ============================================================================
print(f"\n{'='*70}")
print("Analysis Complete!")
print(f"{'='*70}")
print(f"\n📁 Output files:")
print("  - rmsf_data.txt, rmsf_plot.png")
print("  - hbond_data.txt, hbond_plot.png")
print("\n✅ RMSF and H-bond analyses completed successfully!")
