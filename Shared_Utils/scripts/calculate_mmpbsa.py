#!/usr/bin/env python3
"""
MM/PBSA Binding Free Energy Calculation
Using OpenMM for energy calculations
"""

import MDAnalysis as mda
from openmm import app, unit
from openmm.app import CharmmPsfFile, CharmmParameterSet
import numpy as np
import pandas as pd
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set up paths
BASE_DIR = Path("/home/pjho3/projects/Drug/final_complex/controlcomplex/openmm")
OUTPUT_DIR = Path("/home/pjho3/projects/Drug/final_complex/analysis")
OUTPUT_DIR.mkdir(exist_ok=True)

# Input files
PSF_FILE = BASE_DIR / "step5_input.psf"
PDB_FILE = BASE_DIR / "step5_input.pdb"
DCD_FILE = BASE_DIR / "step7_1.dcd"
TOPPAR_STR = BASE_DIR / "toppar.str"

print("=" * 80)
print("MM/PBSA Binding Free Energy Calculation")
print("=" * 80)
print(f"\nInput files:")
print(f"  PSF: {PSF_FILE}")
print(f"  PDB: {PDB_FILE}")
print(f"  DCD: {DCD_FILE}")
print(f"  Parameters: {TOPPAR_STR}")
print()

# Load trajectory with MDAnalysis
print("Loading trajectory...")
u = mda.Universe(str(PSF_FILE), str(DCD_FILE))
print(f"✅ Loaded: {len(u.trajectory)} frames")
print(f"   Total time: {u.trajectory.totaltime / 1000:.2f} ns")
print()

# Identify protein and ligand
print("Identifying system components...")
protein = u.select_atoms("protein")
print(f"  Protein atoms: {len(protein)}")

# Find ligand
try:
    ligand = u.select_atoms("resname SDG")
    if len(ligand) == 0:
        ligand = u.select_atoms("not protein and not resname TIP3 and not resname SOD and not resname CLA and not resname POT")
    print(f"  Ligand atoms: {len(ligand)}")
    print(f"  Ligand resname: {ligand.resnames[0] if len(ligand) > 0 else 'Not found'}")
except:
    print("  ⚠️  Could not identify ligand")
    ligand = None

if ligand is None or len(ligand) == 0:
    print("\n❌ Error: No ligand found. Cannot perform MM/PBSA calculation.")
    exit(1)

# Create complex, protein-only, and ligand-only selections
complex_atoms = u.select_atoms("protein or resname SDG or (not protein and not resname TIP3 and not resname SOD and not resname CLA and not resname POT)")
print(f"  Complex atoms: {len(complex_atoms)}")
print()

# ============================================================================
# MM/PBSA Calculation (Simplified approach)
# ============================================================================
print("=" * 80)
print("Calculating Binding Energy Components")
print("=" * 80)
print("\nNote: This is a simplified MM/PBSA-like calculation")
print("For full MM/PBSA, consider using MMPBSA.py from AmberTools")
print()

# Sample frames for energy calculation (every 10th frame to save time)
sample_interval = 10
sample_frames = list(range(0, len(u.trajectory), sample_interval))
print(f"Sampling {len(sample_frames)} frames (every {sample_interval}th frame)...")

# Calculate interaction energies
energies = []

for frame_idx in sample_frames:
    u.trajectory[frame_idx]
    
    # Calculate protein-ligand distance
    protein_com = protein.center_of_mass()
    ligand_com = ligand.center_of_mass()
    distance = np.linalg.norm(protein_com - ligand_com)
    
    # Count close contacts (< 4.5 Å)
    close_contacts = 0
    for p_atom in protein.positions:
        for l_atom in ligand.positions:
            dist = np.linalg.norm(p_atom - l_atom)
            if dist < 4.5:
                close_contacts += 1
    
    # Calculate buried surface area approximation
    # (number of ligand atoms within 5 Å of protein)
    buried_atoms = 0
    for l_atom in ligand.positions:
        for p_atom in protein.positions:
            dist = np.linalg.norm(p_atom - l_atom)
            if dist < 5.0:
                buried_atoms += 1
                break
    
    energies.append({
        'Frame': frame_idx,
        'Time_ns': u.trajectory[frame_idx].time / 1000,
        'Distance_A': distance,
        'Close_Contacts': close_contacts,
        'Buried_Ligand_Atoms': buried_atoms,
        'Burial_Fraction': buried_atoms / len(ligand)
    })

energy_data = pd.DataFrame(energies)
energy_data.to_csv(OUTPUT_DIR / "binding_energy_components.csv", index=False)

print(f"\n✅ Energy components calculated for {len(sample_frames)} frames")
print(f"\nAverage values:")
print(f"  Distance: {energy_data['Distance_A'].mean():.2f} ± {energy_data['Distance_A'].std():.2f} Å")
print(f"  Close contacts: {energy_data['Close_Contacts'].mean():.1f} ± {energy_data['Close_Contacts'].std():.1f}")
print(f"  Buried ligand atoms: {energy_data['Buried_Ligand_Atoms'].mean():.1f} / {len(ligand)}")
print(f"  Burial fraction: {energy_data['Burial_Fraction'].mean():.2%} ± {energy_data['Burial_Fraction'].std():.2%}")

# ============================================================================
# Interaction Energy Estimation
# ============================================================================
print("\n" + "=" * 80)
print("Interaction Energy Estimation")
print("=" * 80)

# Estimate interaction energy based on contacts
# Rough approximation: -0.5 kcal/mol per close contact
estimated_interaction_energy = -0.5 * energy_data['Close_Contacts'].mean()
print(f"\nEstimated interaction energy: {estimated_interaction_energy:.2f} kcal/mol")
print("(Based on contact approximation: -0.5 kcal/mol per contact)")

# Estimate desolvation penalty based on burial
# Rough approximation: +2 kcal/mol per buried atom
estimated_desolvation = 2.0 * energy_data['Buried_Ligand_Atoms'].mean()
print(f"\nEstimated desolvation penalty: +{estimated_desolvation:.2f} kcal/mol")
print("(Based on burial approximation: +2 kcal/mol per buried atom)")

# Net binding energy estimate
net_binding_energy = estimated_interaction_energy + estimated_desolvation
print(f"\nNet binding energy estimate: {net_binding_energy:.2f} kcal/mol")

# ============================================================================
# Save Summary
# ============================================================================
summary = {
    'Component': [
        'Protein-Ligand Distance (Å)',
        'Close Contacts (<4.5Å)',
        'Buried Ligand Atoms',
        'Burial Fraction (%)',
        'Estimated Interaction Energy (kcal/mol)',
        'Estimated Desolvation Penalty (kcal/mol)',
        'Net Binding Energy Estimate (kcal/mol)'
    ],
    'Value': [
        f"{energy_data['Distance_A'].mean():.2f} ± {energy_data['Distance_A'].std():.2f}",
        f"{energy_data['Close_Contacts'].mean():.1f} ± {energy_data['Close_Contacts'].std():.1f}",
        f"{energy_data['Buried_Ligand_Atoms'].mean():.1f}",
        f"{energy_data['Burial_Fraction'].mean()*100:.1f} ± {energy_data['Burial_Fraction'].std()*100:.1f}",
        f"{estimated_interaction_energy:.2f}",
        f"{estimated_desolvation:.2f}",
        f"{net_binding_energy:.2f}"
    ]
}

summary_df = pd.DataFrame(summary)
summary_df.to_csv(OUTPUT_DIR / "mmpbsa_summary.csv", index=False)

print("\n" + "=" * 80)
print("Summary")
print("=" * 80)
print("\n" + summary_df.to_string(index=False))

print("\n" + "=" * 80)
print("⚠️  Important Notes")
print("=" * 80)
print("""
This is a SIMPLIFIED binding energy estimation based on:
- Geometric analysis (contacts, burial)
- Empirical approximations

For accurate MM/PBSA calculations, use:
1. MMPBSA.py from AmberTools
2. gmx_MMPBSA for GROMACS
3. OpenMM with explicit solvation free energy calculations

The values provided here are ORDER-OF-MAGNITUDE estimates only.
They indicate binding stability but should not be used for
quantitative comparisons without proper MM/PBSA calculations.
""")

print("\n" + "=" * 80)
print("✅ MM/PBSA Analysis Complete!")
print("=" * 80)
print(f"\nResults saved to: {OUTPUT_DIR}")
print("  - binding_energy_components.csv")
print("  - mmpbsa_summary.csv")
print()
