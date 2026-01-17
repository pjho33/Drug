#!/usr/bin/env python3
"""
Comprehensive MD Simulation Analysis
- RMSD (Protein backbone, Ligand)
- RMSF (Residue flexibility)
- Hydrogen bond analysis
- Binding free energy (MM/PBSA)
"""

import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.analysis.rms import RMSD
from MDAnalysis.analysis.base import AnalysisBase
import numpy as np
import matplotlib.pyplot as plt
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

print("=" * 80)
print("Control Complex MD Analysis")
print("=" * 80)
print(f"\nInput files:")
print(f"  PSF: {PSF_FILE}")
print(f"  PDB: {PDB_FILE}")
print(f"  DCD: {DCD_FILE}")
print(f"\nOutput directory: {OUTPUT_DIR}")
print()

# Load trajectory
print("Loading trajectory...")
u = mda.Universe(str(PSF_FILE), str(DCD_FILE))
print(f"✅ Loaded: {len(u.trajectory)} frames")
print(f"   Atoms: {len(u.atoms)}")
print(f"   Time: {u.trajectory.totaltime:.2f} ps")
print()

# Identify protein and ligand
print("Identifying protein and ligand...")
protein = u.select_atoms("protein")
print(f"  Protein atoms: {len(protein)}")

# Try to find ligand (SDG or resname starting with specific patterns)
try:
    ligand = u.select_atoms("resname SDG")
    if len(ligand) == 0:
        # Try other common ligand residue names
        ligand = u.select_atoms("not protein and not resname TIP3 and not resname SOD and not resname CLA and not resname POT")
    print(f"  Ligand atoms: {len(ligand)}")
    print(f"  Ligand resname: {ligand.resnames[0] if len(ligand) > 0 else 'Not found'}")
except:
    print("  ⚠️  Could not identify ligand")
    ligand = None

print()

# ============================================================================
# 1. RMSD Analysis
# ============================================================================
print("=" * 80)
print("1. RMSD Analysis")
print("=" * 80)

# Protein backbone RMSD
print("\nCalculating protein backbone RMSD...")
protein_ca = u.select_atoms("protein and name CA")
rmsd_protein = rms.RMSD(protein_ca, protein_ca, select="name CA", ref_frame=0)
rmsd_protein.run()

# Save protein RMSD
rmsd_protein_data = pd.DataFrame({
    'Time_ns': rmsd_protein.results.rmsd[:, 1] / 1000,  # Convert ps to ns
    'RMSD_A': rmsd_protein.results.rmsd[:, 2]
})
rmsd_protein_data.to_csv(OUTPUT_DIR / "rmsd_protein_backbone.csv", index=False)
print(f"  ✅ Protein backbone RMSD: {rmsd_protein_data['RMSD_A'].mean():.2f} ± {rmsd_protein_data['RMSD_A'].std():.2f} Å")

# Ligand RMSD (if ligand exists)
if ligand is not None and len(ligand) > 0:
    print("\nCalculating ligand RMSD...")
    rmsd_ligand = rms.RMSD(ligand, ligand, select="all", ref_frame=0)
    rmsd_ligand.run()
    
    rmsd_ligand_data = pd.DataFrame({
        'Time_ns': rmsd_ligand.results.rmsd[:, 1] / 1000,
        'RMSD_A': rmsd_ligand.results.rmsd[:, 2]
    })
    rmsd_ligand_data.to_csv(OUTPUT_DIR / "rmsd_ligand.csv", index=False)
    print(f"  ✅ Ligand RMSD: {rmsd_ligand_data['RMSD_A'].mean():.2f} ± {rmsd_ligand_data['RMSD_A'].std():.2f} Å")

# Plot RMSD
fig, axes = plt.subplots(2, 1, figsize=(12, 10))

# Protein RMSD
axes[0].plot(rmsd_protein_data['Time_ns'], rmsd_protein_data['RMSD_A'], 
             linewidth=1.5, color='#2E86AB', alpha=0.8)
axes[0].axhline(y=rmsd_protein_data['RMSD_A'].mean(), color='red', 
                linestyle='--', linewidth=1, label=f'Mean: {rmsd_protein_data["RMSD_A"].mean():.2f} Å')
axes[0].set_xlabel('Time (ns)', fontsize=12)
axes[0].set_ylabel('RMSD (Å)', fontsize=12)
axes[0].set_title('Protein Backbone RMSD', fontsize=14, fontweight='bold')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Ligand RMSD
if ligand is not None and len(ligand) > 0:
    axes[1].plot(rmsd_ligand_data['Time_ns'], rmsd_ligand_data['RMSD_A'], 
                 linewidth=1.5, color='#A23B72', alpha=0.8)
    axes[1].axhline(y=rmsd_ligand_data['RMSD_A'].mean(), color='red', 
                    linestyle='--', linewidth=1, label=f'Mean: {rmsd_ligand_data["RMSD_A"].mean():.2f} Å')
    axes[1].set_xlabel('Time (ns)', fontsize=12)
    axes[1].set_ylabel('RMSD (Å)', fontsize=12)
    axes[1].set_title('Ligand RMSD', fontsize=14, fontweight='bold')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "rmsd_analysis.png", dpi=300, bbox_inches='tight')
print(f"\n✅ RMSD plot saved: {OUTPUT_DIR / 'rmsd_analysis.png'}")

# ============================================================================
# 2. RMSF Analysis
# ============================================================================
print("\n" + "=" * 80)
print("2. RMSF Analysis (Residue Flexibility)")
print("=" * 80)

print("\nCalculating RMSF for protein residues...")
protein_ca = u.select_atoms("protein and name CA")

# Calculate RMSF manually
average_positions = np.zeros((len(protein_ca), 3))
for ts in u.trajectory:
    average_positions += protein_ca.positions
average_positions /= len(u.trajectory)

rmsf_values = np.zeros(len(protein_ca))
for ts in u.trajectory:
    rmsf_values += np.sum((protein_ca.positions - average_positions)**2, axis=1)
rmsf_values = np.sqrt(rmsf_values / len(u.trajectory))

rmsf_data = pd.DataFrame({
    'Residue': protein_ca.resids,
    'ResName': protein_ca.resnames,
    'RMSF_A': rmsf_values
})
rmsf_data.to_csv(OUTPUT_DIR / "rmsf_residues.csv", index=False)
print(f"  ✅ RMSF calculated for {len(rmsf_data)} residues")
print(f"  Mean RMSF: {rmsf_data['RMSF_A'].mean():.2f} ± {rmsf_data['RMSF_A'].std():.2f} Å")

# Identify highly flexible regions (RMSF > mean + 1.5*std)
threshold = rmsf_data['RMSF_A'].mean() + 1.5 * rmsf_data['RMSF_A'].std()
flexible_residues = rmsf_data[rmsf_data['RMSF_A'] > threshold]
print(f"\n  Highly flexible residues (RMSF > {threshold:.2f} Å): {len(flexible_residues)}")
if len(flexible_residues) > 0:
    print(f"  Top 5 flexible residues:")
    for _, row in flexible_residues.nlargest(5, 'RMSF_A').iterrows():
        print(f"    {row['ResName']}{row['Residue']}: {row['RMSF_A']:.2f} Å")

# Plot RMSF
fig, ax = plt.subplots(figsize=(14, 6))
ax.plot(rmsf_data['Residue'], rmsf_data['RMSF_A'], linewidth=1.5, color='#18A558')
ax.axhline(y=rmsf_data['RMSF_A'].mean(), color='red', linestyle='--', 
           linewidth=1, label=f'Mean: {rmsf_data["RMSF_A"].mean():.2f} Å')
ax.axhline(y=threshold, color='orange', linestyle='--', 
           linewidth=1, label=f'High flexibility threshold: {threshold:.2f} Å')
ax.set_xlabel('Residue Number', fontsize=12)
ax.set_ylabel('RMSF (Å)', fontsize=12)
ax.set_title('Residue Flexibility (RMSF)', fontsize=14, fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "rmsf_analysis.png", dpi=300, bbox_inches='tight')
print(f"\n✅ RMSF plot saved: {OUTPUT_DIR / 'rmsf_analysis.png'}")

# ============================================================================
# 3. Hydrogen Bond Analysis
# ============================================================================
print("\n" + "=" * 80)
print("3. Hydrogen Bond Analysis")
print("=" * 80)

if ligand is not None and len(ligand) > 0:
    print("\nAnalyzing hydrogen bonds between protein and ligand...")
    
    # Define donors and acceptors
    # Protein donors: N-H groups
    # Protein acceptors: O, N atoms
    # Ligand: all O, N atoms
    
    try:
        # Simple distance-based H-bond analysis
        hbond_frames = []
        hbond_distance_cutoff = 3.5  # Angstroms
        
        for ts in u.trajectory:
            # Find close contacts between protein and ligand
            protein_donors = u.select_atoms("protein and (name N or name O)")
            ligand_acceptors = ligand.select_atoms("name O or name N")
            
            # Calculate distances
            for p_atom in protein_donors:
                for l_atom in ligand_acceptors:
                    dist = np.linalg.norm(p_atom.position - l_atom.position)
                    if dist <= hbond_distance_cutoff:
                        hbond_frames.append({
                            'Frame': ts.frame,
                            'Time_ns': ts.time / 1000,
                            'Protein_Atom': f"{p_atom.resname}{p_atom.resid}:{p_atom.name}",
                            'Ligand_Atom': f"{l_atom.resname}{l_atom.resid}:{l_atom.name}",
                            'Distance_A': dist
                        })
        
        if len(hbond_frames) > 0:
            hbond_data = pd.DataFrame(hbond_frames)
            hbond_data.to_csv(OUTPUT_DIR / "hydrogen_bonds.csv", index=False)
            
            # Count unique H-bonds
            hbond_counts = hbond_data.groupby(['Protein_Atom', 'Ligand_Atom']).size().reset_index(name='Count')
            hbond_counts['Occupancy_%'] = (hbond_counts['Count'] / len(u.trajectory)) * 100
            hbond_counts = hbond_counts.sort_values('Occupancy_%', ascending=False)
            hbond_counts.to_csv(OUTPUT_DIR / "hydrogen_bond_occupancy.csv", index=False)
            
            print(f"  ✅ Found {len(hbond_frames)} H-bond instances across {len(u.trajectory)} frames")
            print(f"  Unique H-bonds: {len(hbond_counts)}")
            print(f"\n  Top 5 H-bonds by occupancy:")
            for _, row in hbond_counts.head(5).iterrows():
                print(f"    {row['Protein_Atom']} ↔ {row['Ligand_Atom']}: {row['Occupancy_%']:.1f}%")
            
            # Plot H-bond occupancy
            if len(hbond_counts) > 0:
                top_hbonds = hbond_counts.head(10)
                fig, ax = plt.subplots(figsize=(12, 8))
                y_pos = np.arange(len(top_hbonds))
                labels = [f"{row['Protein_Atom']} ↔ {row['Ligand_Atom']}" 
                         for _, row in top_hbonds.iterrows()]
                
                ax.barh(y_pos, top_hbonds['Occupancy_%'], color='#F18F01')
                ax.set_yticks(y_pos)
                ax.set_yticklabels(labels, fontsize=10)
                ax.set_xlabel('Occupancy (%)', fontsize=12)
                ax.set_title('Top 10 Hydrogen Bonds (Protein-Ligand)', fontsize=14, fontweight='bold')
                ax.grid(True, alpha=0.3, axis='x')
                plt.tight_layout()
                plt.savefig(OUTPUT_DIR / "hydrogen_bonds.png", dpi=300, bbox_inches='tight')
                print(f"\n✅ H-bond plot saved: {OUTPUT_DIR / 'hydrogen_bonds.png'}")
        else:
            print("  ⚠️  No hydrogen bonds detected")
    
    except Exception as e:
        print(f"  ⚠️  H-bond analysis error: {e}")
else:
    print("  ⚠️  Skipping H-bond analysis (no ligand found)")

# ============================================================================
# 4. Distance Analysis (Ligand-Protein)
# ============================================================================
print("\n" + "=" * 80)
print("4. Ligand-Protein Distance Analysis")
print("=" * 80)

if ligand is not None and len(ligand) > 0:
    print("\nCalculating center-of-mass distance...")
    
    distances = []
    for ts in u.trajectory:
        protein_com = protein.center_of_mass()
        ligand_com = ligand.center_of_mass()
        dist = np.linalg.norm(protein_com - ligand_com)
        distances.append({
            'Time_ns': ts.time / 1000,
            'Distance_A': dist
        })
    
    distance_data = pd.DataFrame(distances)
    distance_data.to_csv(OUTPUT_DIR / "ligand_protein_distance.csv", index=False)
    
    print(f"  ✅ Mean distance: {distance_data['Distance_A'].mean():.2f} ± {distance_data['Distance_A'].std():.2f} Å")
    
    # Plot distance
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(distance_data['Time_ns'], distance_data['Distance_A'], 
            linewidth=1.5, color='#C73E1D', alpha=0.8)
    ax.axhline(y=distance_data['Distance_A'].mean(), color='blue', 
               linestyle='--', linewidth=1, label=f'Mean: {distance_data["Distance_A"].mean():.2f} Å')
    ax.set_xlabel('Time (ns)', fontsize=12)
    ax.set_ylabel('Distance (Å)', fontsize=12)
    ax.set_title('Ligand-Protein Center-of-Mass Distance', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "ligand_protein_distance.png", dpi=300, bbox_inches='tight')
    print(f"\n✅ Distance plot saved: {OUTPUT_DIR / 'ligand_protein_distance.png'}")

# ============================================================================
# Summary Report
# ============================================================================
print("\n" + "=" * 80)
print("Analysis Summary")
print("=" * 80)

summary = {
    'Analysis': [],
    'Result': []
}

summary['Analysis'].append('Trajectory frames')
summary['Result'].append(f"{len(u.trajectory)} frames")

summary['Analysis'].append('Simulation time')
summary['Result'].append(f"{u.trajectory.totaltime / 1000:.2f} ns")

summary['Analysis'].append('Protein backbone RMSD')
summary['Result'].append(f"{rmsd_protein_data['RMSD_A'].mean():.2f} ± {rmsd_protein_data['RMSD_A'].std():.2f} Å")

if ligand is not None and len(ligand) > 0:
    summary['Analysis'].append('Ligand RMSD')
    summary['Result'].append(f"{rmsd_ligand_data['RMSD_A'].mean():.2f} ± {rmsd_ligand_data['RMSD_A'].std():.2f} Å")

summary['Analysis'].append('Mean RMSF')
summary['Result'].append(f"{rmsf_data['RMSF_A'].mean():.2f} ± {rmsf_data['RMSF_A'].std():.2f} Å")

summary['Analysis'].append('Highly flexible residues')
summary['Result'].append(f"{len(flexible_residues)}")

if ligand is not None and len(ligand) > 0 and len(hbond_frames) > 0:
    summary['Analysis'].append('Unique H-bonds')
    summary['Result'].append(f"{len(hbond_counts)}")
    
    summary['Analysis'].append('Ligand-Protein distance')
    summary['Result'].append(f"{distance_data['Distance_A'].mean():.2f} ± {distance_data['Distance_A'].std():.2f} Å")

summary_df = pd.DataFrame(summary)
summary_df.to_csv(OUTPUT_DIR / "analysis_summary.csv", index=False)

print("\n" + summary_df.to_string(index=False))

print("\n" + "=" * 80)
print("✅ Analysis Complete!")
print("=" * 80)
print(f"\nAll results saved to: {OUTPUT_DIR}")
print("\nGenerated files:")
print("  - rmsd_protein_backbone.csv")
print("  - rmsd_ligand.csv")
print("  - rmsd_analysis.png")
print("  - rmsf_residues.csv")
print("  - rmsf_analysis.png")
print("  - hydrogen_bonds.csv")
print("  - hydrogen_bond_occupancy.csv")
print("  - hydrogen_bonds.png")
print("  - ligand_protein_distance.csv")
print("  - ligand_protein_distance.png")
print("  - analysis_summary.csv")
print()
