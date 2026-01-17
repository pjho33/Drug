#!/usr/bin/env python
"""
MD Simulation Analysis
======================
Analyzes RMSD, energy, and ligand-receptor interactions
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import pandas as pd

# Paths
RESULTS_DIR = "/home/pjho3/projects/Drug/results/phase2_gnina_md/model_b"
RECEPTOR_PDB = "/home/pjho3/projects/Drug/results/phase2_gnina_md/cleaned_receptor.pdb"
OUTPUT_DIR = "/home/pjho3/projects/Drug/results/analysis"

os.makedirs(OUTPUT_DIR, exist_ok=True)


def analyze_energy_log(log_path):
    """Analyze energy and temperature from simulation log"""
    print("\n" + "=" * 60)
    print("Energy Analysis")
    print("=" * 60)
    
    df = pd.read_csv(log_path)
    df.columns = ['Step', 'Potential_Energy', 'Temperature', 'Volume', 'Speed']
    
    # Convert to ns
    dt_fs = 2.0  # femtoseconds per step
    df['Time_ns'] = df['Step'] * dt_fs / 1e6
    
    # Statistics
    print(f"\nSimulation Statistics:")
    print(f"  Total time: {df['Time_ns'].max():.2f} ns")
    print(f"  Total steps: {df['Step'].max()}")
    print(f"\nTemperature:")
    print(f"  Mean: {df['Temperature'].mean():.2f} K")
    print(f"  Std: {df['Temperature'].std():.2f} K")
    print(f"\nPotential Energy:")
    print(f"  Mean: {df['Potential_Energy'].mean()/1000:.2f} MJ/mol")
    print(f"  Std: {df['Potential_Energy'].std()/1000:.2f} MJ/mol")
    
    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    
    # Temperature
    axes[0, 0].plot(df['Time_ns'], df['Temperature'], 'b-', alpha=0.7)
    axes[0, 0].axhline(y=310, color='r', linestyle='--', label='Target (310K)')
    axes[0, 0].set_xlabel('Time (ns)')
    axes[0, 0].set_ylabel('Temperature (K)')
    axes[0, 0].set_title('Temperature vs Time')
    axes[0, 0].legend()
    
    # Potential Energy
    axes[0, 1].plot(df['Time_ns'], df['Potential_Energy']/1000, 'g-', alpha=0.7)
    axes[0, 1].set_xlabel('Time (ns)')
    axes[0, 1].set_ylabel('Potential Energy (MJ/mol)')
    axes[0, 1].set_title('Potential Energy vs Time')
    
    # Volume
    axes[1, 0].plot(df['Time_ns'], df['Volume'], 'm-', alpha=0.7)
    axes[1, 0].set_xlabel('Time (ns)')
    axes[1, 0].set_ylabel('Box Volume (nm³)')
    axes[1, 0].set_title('Box Volume vs Time')
    
    # Speed
    axes[1, 1].plot(df['Time_ns'], df['Speed'], 'c-', alpha=0.7)
    axes[1, 1].set_xlabel('Time (ns)')
    axes[1, 1].set_ylabel('Speed (ns/day)')
    axes[1, 1].set_title('Simulation Speed')
    
    plt.tight_layout()
    plot_path = os.path.join(OUTPUT_DIR, 'energy_analysis.png')
    plt.savefig(plot_path, dpi=150)
    print(f"\n✅ Energy plot saved: {plot_path}")
    plt.close()
    
    return df


def analyze_rmsd(traj_path, topology_path):
    """Calculate RMSD for protein and ligand"""
    print("\n" + "=" * 60)
    print("RMSD Analysis")
    print("=" * 60)
    
    print(f"Loading trajectory: {traj_path}")
    print(f"This may take a while for large trajectories...")
    
    # Load trajectory (sample every 10 frames to reduce memory)
    traj = md.load(traj_path, top=topology_path, stride=10)
    print(f"Loaded {traj.n_frames} frames, {traj.n_atoms} atoms")
    
    # Get time in ns
    dt_fs = 2.0
    save_ps = 10.0
    stride = 10
    time_ns = np.arange(traj.n_frames) * (save_ps / 1000) * stride
    
    # Select protein backbone (CA atoms)
    protein_atoms = traj.topology.select('protein and name CA')
    print(f"Protein CA atoms: {len(protein_atoms)}")
    
    if len(protein_atoms) > 0:
        # Calculate protein RMSD
        protein_rmsd = md.rmsd(traj, traj, frame=0, atom_indices=protein_atoms) * 10  # nm to Å
        
        print(f"\nProtein Backbone RMSD:")
        print(f"  Initial: {protein_rmsd[0]:.2f} Å")
        print(f"  Final: {protein_rmsd[-1]:.2f} Å")
        print(f"  Mean: {protein_rmsd.mean():.2f} Å")
        print(f"  Max: {protein_rmsd.max():.2f} Å")
    else:
        protein_rmsd = None
        print("⚠️ No protein CA atoms found")
    
    # Select ligand (non-protein, non-water)
    ligand_atoms = traj.topology.select('not protein and not water and not (name Na or name Cl)')
    print(f"Ligand atoms: {len(ligand_atoms)}")
    
    if len(ligand_atoms) > 0:
        # Calculate ligand RMSD (after aligning on protein)
        if len(protein_atoms) > 0:
            # Superpose on protein first
            traj_aligned = traj.superpose(traj, frame=0, atom_indices=protein_atoms)
            ligand_rmsd = md.rmsd(traj_aligned, traj_aligned, frame=0, atom_indices=ligand_atoms) * 10
        else:
            ligand_rmsd = md.rmsd(traj, traj, frame=0, atom_indices=ligand_atoms) * 10
        
        print(f"\nLigand RMSD (aligned on protein):")
        print(f"  Initial: {ligand_rmsd[0]:.2f} Å")
        print(f"  Final: {ligand_rmsd[-1]:.2f} Å")
        print(f"  Mean: {ligand_rmsd.mean():.2f} Å")
        print(f"  Max: {ligand_rmsd.max():.2f} Å")
    else:
        ligand_rmsd = None
        print("⚠️ No ligand atoms found")
    
    # Plot RMSD
    fig, ax = plt.subplots(figsize=(10, 6))
    
    if protein_rmsd is not None:
        ax.plot(time_ns, protein_rmsd, 'b-', label='Protein Backbone', alpha=0.8)
    if ligand_rmsd is not None:
        ax.plot(time_ns, ligand_rmsd, 'r-', label='Ligand', alpha=0.8)
    
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('RMSD (Å)')
    ax.set_title('RMSD vs Time - Model B (Arg-TRIS-PEG2)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plot_path = os.path.join(OUTPUT_DIR, 'rmsd_analysis.png')
    plt.savefig(plot_path, dpi=150)
    print(f"\n✅ RMSD plot saved: {plot_path}")
    plt.close()
    
    return protein_rmsd, ligand_rmsd, time_ns


def analyze_contacts(traj_path, topology_path):
    """Analyze protein-ligand contacts"""
    print("\n" + "=" * 60)
    print("Contact Analysis")
    print("=" * 60)
    
    # Load trajectory (sample every 10 frames)
    traj = md.load(traj_path, top=topology_path, stride=10)
    
    # Get ligand and protein atoms
    ligand_atoms = traj.topology.select('not protein and not water and not (name Na or name Cl)')
    protein_atoms = traj.topology.select('protein')
    
    if len(ligand_atoms) == 0 or len(protein_atoms) == 0:
        print("⚠️ Cannot analyze contacts - missing atoms")
        return None
    
    # Calculate contacts (within 4 Å)
    cutoff = 0.4  # nm
    contacts_per_frame = []
    
    for i in range(traj.n_frames):
        frame = traj[i]
        ligand_pos = frame.xyz[0, ligand_atoms]
        protein_pos = frame.xyz[0, protein_atoms]
        
        # Calculate distances
        n_contacts = 0
        for lig_pos in ligand_pos:
            dists = np.sqrt(np.sum((protein_pos - lig_pos)**2, axis=1))
            n_contacts += np.sum(dists < cutoff)
        
        contacts_per_frame.append(n_contacts)
    
    contacts = np.array(contacts_per_frame)
    
    print(f"\nProtein-Ligand Contacts (< 4 Å):")
    print(f"  Initial: {contacts[0]}")
    print(f"  Final: {contacts[-1]}")
    print(f"  Mean: {contacts.mean():.1f}")
    print(f"  Min: {contacts.min()}")
    print(f"  Max: {contacts.max()}")
    
    # Plot
    dt_fs = 2.0
    save_ps = 10.0
    stride = 10
    time_ns = np.arange(len(contacts)) * (save_ps / 1000) * stride
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(time_ns, contacts, 'g-', alpha=0.8)
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Number of Contacts')
    ax.set_title('Protein-Ligand Contacts (< 4 Å) - Model B')
    ax.grid(True, alpha=0.3)
    
    plot_path = os.path.join(OUTPUT_DIR, 'contacts_analysis.png')
    plt.savefig(plot_path, dpi=150)
    print(f"\n✅ Contacts plot saved: {plot_path}")
    plt.close()
    
    return contacts


def main():
    print("=" * 60)
    print("MD Simulation Analysis - Model B (Arg-TRIS-PEG2)")
    print("=" * 60)
    
    # File paths
    log_path = os.path.join(RESULTS_DIR, "prod_model_b_log.csv")
    traj_path = os.path.join(RESULTS_DIR, "prod_model_b.dcd")
    final_pdb = os.path.join(RESULTS_DIR, "prod_model_b_final.pdb")
    
    # 1. Energy analysis
    energy_df = analyze_energy_log(log_path)
    
    # 2. RMSD analysis
    protein_rmsd, ligand_rmsd, time_ns = analyze_rmsd(traj_path, final_pdb)
    
    # 3. Contact analysis
    contacts = analyze_contacts(traj_path, final_pdb)
    
    # Summary
    print("\n" + "=" * 60)
    print("Analysis Summary")
    print("=" * 60)
    print(f"\nSimulation: Model B (Arg-TRIS-PEG2) with GLUT1")
    print(f"Duration: {energy_df['Time_ns'].max():.2f} ns")
    print(f"\nStability Metrics:")
    print(f"  Temperature: {energy_df['Temperature'].mean():.1f} ± {energy_df['Temperature'].std():.1f} K")
    if protein_rmsd is not None:
        print(f"  Protein RMSD: {protein_rmsd.mean():.2f} ± {protein_rmsd.std():.2f} Å")
    if ligand_rmsd is not None:
        print(f"  Ligand RMSD: {ligand_rmsd.mean():.2f} ± {ligand_rmsd.std():.2f} Å")
    if contacts is not None:
        print(f"  Contacts: {contacts.mean():.0f} ± {contacts.std():.0f}")
    
    print(f"\nOutput directory: {OUTPUT_DIR}")
    print("✅ Analysis complete!")


if __name__ == "__main__":
    main()
