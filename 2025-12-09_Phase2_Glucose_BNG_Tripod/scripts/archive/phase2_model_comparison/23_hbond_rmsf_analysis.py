#!/usr/bin/env python
"""
Hydrogen Bond and RMSF Analysis
================================
Analyzes H-bonds between protein and ligand, and residue fluctuations
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import json
import warnings
warnings.filterwarnings('ignore')

# Paths
RESULTS_DIR = "/home/pjho3/projects/Drug/results/phase2_gnina_md/model_b"
OUTPUT_DIR = "/home/pjho3/projects/Drug/results/analysis"

os.makedirs(OUTPUT_DIR, exist_ok=True)


def analyze_hydrogen_bonds(traj_path, topology_path):
    """
    Analyze hydrogen bonds between protein and ligand
    """
    print("=" * 60)
    print("Hydrogen Bond Analysis")
    print("=" * 60)
    
    # Load trajectory
    print("\nLoading trajectory...")
    traj = md.load(traj_path, top=topology_path, stride=10)
    print(f"Loaded {traj.n_frames} frames, {traj.n_atoms} atoms")
    
    # Get protein and ligand atom indices
    protein_idx = traj.topology.select('protein')
    ligand_idx = traj.topology.select('not protein and not water and not (name Na or name Cl)')
    
    print(f"Protein atoms: {len(protein_idx)}")
    print(f"Ligand atoms: {len(ligand_idx)}")
    
    # Find hydrogen bonds using MDTraj's baker_hubbard method
    # This finds D-H...A where D and A are N, O, or S
    print("\nFinding hydrogen bonds...")
    
    # Get all H-bonds in the system
    hbonds_all = md.baker_hubbard(traj, freq=0.1)  # H-bonds present in at least 10% of frames
    
    # Filter for protein-ligand H-bonds
    protein_set = set(protein_idx)
    ligand_set = set(ligand_idx)
    
    protein_ligand_hbonds = []
    for hbond in hbonds_all:
        donor, hydrogen, acceptor = hbond
        
        # Check if one is in protein and one is in ligand
        donor_in_protein = donor in protein_set
        donor_in_ligand = donor in ligand_set
        acceptor_in_protein = acceptor in protein_set
        acceptor_in_ligand = acceptor in ligand_set
        
        if (donor_in_protein and acceptor_in_ligand) or (donor_in_ligand and acceptor_in_protein):
            protein_ligand_hbonds.append(hbond)
    
    print(f"\nProtein-Ligand H-bonds found: {len(protein_ligand_hbonds)}")
    
    # Analyze each H-bond
    hbond_details = []
    if len(protein_ligand_hbonds) > 0:
        print("\nH-bond details:")
        print("-" * 70)
        
        for i, hbond in enumerate(protein_ligand_hbonds):
            donor, hydrogen, acceptor = hbond
            
            # Get atom info
            donor_atom = traj.topology.atom(donor)
            acceptor_atom = traj.topology.atom(acceptor)
            
            # Calculate occupancy (fraction of frames with this H-bond)
            distances = md.compute_distances(traj, [[hydrogen, acceptor]])
            angles = md.compute_angles(traj, [[donor, hydrogen, acceptor]])
            
            # H-bond criteria: distance < 2.5 Å, angle > 120°
            hbond_present = (distances[:, 0] < 0.25) & (angles[:, 0] > np.radians(120))
            occupancy = np.mean(hbond_present) * 100
            
            mean_dist = np.mean(distances[hbond_present, 0]) * 10 if np.any(hbond_present) else 0
            
            detail = {
                'donor': f"{donor_atom.residue.name}{donor_atom.residue.resSeq}-{donor_atom.name}",
                'acceptor': f"{acceptor_atom.residue.name}{acceptor_atom.residue.resSeq}-{acceptor_atom.name}",
                'occupancy': occupancy,
                'mean_distance': mean_dist
            }
            hbond_details.append(detail)
            
            if occupancy > 5:  # Only show significant H-bonds
                print(f"  {detail['donor']:20s} --> {detail['acceptor']:20s}  "
                      f"Occ: {occupancy:5.1f}%  Dist: {mean_dist:.2f} Å")
    
    # Count H-bonds per frame
    hbond_counts = []
    for frame_idx in range(traj.n_frames):
        frame_hbonds = md.baker_hubbard(traj[frame_idx], freq=0.0)
        
        count = 0
        for hbond in frame_hbonds:
            donor, hydrogen, acceptor = hbond
            donor_in_protein = donor in protein_set
            donor_in_ligand = donor in ligand_set
            acceptor_in_protein = acceptor in protein_set
            acceptor_in_ligand = acceptor in ligand_set
            
            if (donor_in_protein and acceptor_in_ligand) or (donor_in_ligand and acceptor_in_protein):
                count += 1
        
        hbond_counts.append(count)
    
    hbond_counts = np.array(hbond_counts)
    
    print(f"\nH-bond count per frame:")
    print(f"  Mean: {hbond_counts.mean():.1f}")
    print(f"  Std: {hbond_counts.std():.1f}")
    print(f"  Min: {hbond_counts.min()}")
    print(f"  Max: {hbond_counts.max()}")
    
    # Plot H-bond count over time
    dt_fs = 2.0
    save_ps = 10.0
    stride = 10
    time_ns = np.arange(len(hbond_counts)) * (save_ps / 1000) * stride
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(time_ns, hbond_counts, 'b-', alpha=0.8, linewidth=1.5)
    ax.axhline(y=hbond_counts.mean(), color='r', linestyle='--', label=f'Mean: {hbond_counts.mean():.1f}')
    ax.fill_between(time_ns, hbond_counts.mean() - hbond_counts.std(), 
                    hbond_counts.mean() + hbond_counts.std(), alpha=0.2, color='r')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Number of H-bonds')
    ax.set_title('Protein-Ligand Hydrogen Bonds - Model B (Arg-TRIS-PEG2)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plot_path = os.path.join(OUTPUT_DIR, 'hbond_analysis.png')
    plt.savefig(plot_path, dpi=150)
    print(f"\n✅ H-bond plot saved: {plot_path}")
    plt.close()
    
    return hbond_counts, hbond_details


def analyze_rmsf(traj_path, topology_path):
    """
    Calculate Root Mean Square Fluctuation (RMSF) per residue
    """
    print("\n" + "=" * 60)
    print("RMSF Analysis")
    print("=" * 60)
    
    # Load trajectory
    print("\nLoading trajectory...")
    traj = md.load(traj_path, top=topology_path, stride=10)
    print(f"Loaded {traj.n_frames} frames")
    
    # Select protein CA atoms
    ca_atoms = traj.topology.select('protein and name CA')
    print(f"CA atoms: {len(ca_atoms)}")
    
    if len(ca_atoms) == 0:
        print("⚠️ No CA atoms found!")
        return None, None
    
    # Superpose on first frame
    traj.superpose(traj, frame=0, atom_indices=ca_atoms)
    
    # Calculate RMSF
    # RMSF = sqrt(mean((x - mean(x))^2))
    ca_positions = traj.xyz[:, ca_atoms, :]  # (n_frames, n_ca, 3)
    mean_positions = np.mean(ca_positions, axis=0)  # (n_ca, 3)
    
    fluctuations = ca_positions - mean_positions  # (n_frames, n_ca, 3)
    rmsf = np.sqrt(np.mean(np.sum(fluctuations**2, axis=2), axis=0)) * 10  # nm to Å
    
    # Get residue info
    residues = []
    for idx in ca_atoms:
        atom = traj.topology.atom(idx)
        residues.append(f"{atom.residue.name}{atom.residue.resSeq}")
    
    residue_ids = [traj.topology.atom(idx).residue.resSeq for idx in ca_atoms]
    
    print(f"\nRMSF Statistics:")
    print(f"  Mean: {rmsf.mean():.2f} Å")
    print(f"  Std: {rmsf.std():.2f} Å")
    print(f"  Min: {rmsf.min():.2f} Å (Residue {residues[np.argmin(rmsf)]})")
    print(f"  Max: {rmsf.max():.2f} Å (Residue {residues[np.argmax(rmsf)]})")
    
    # Find flexible regions (RMSF > mean + 1*std)
    threshold = rmsf.mean() + rmsf.std()
    flexible_idx = np.where(rmsf > threshold)[0]
    
    print(f"\nFlexible regions (RMSF > {threshold:.2f} Å):")
    if len(flexible_idx) > 0:
        for idx in flexible_idx[:10]:  # Show top 10
            print(f"  {residues[idx]}: {rmsf[idx]:.2f} Å")
    
    # Find rigid regions (RMSF < mean - 0.5*std)
    rigid_threshold = max(0, rmsf.mean() - 0.5*rmsf.std())
    rigid_idx = np.where(rmsf < rigid_threshold)[0]
    
    print(f"\nRigid regions (RMSF < {rigid_threshold:.2f} Å):")
    if len(rigid_idx) > 0:
        for idx in rigid_idx[:10]:  # Show top 10
            print(f"  {residues[idx]}: {rmsf[idx]:.2f} Å")
    
    # Plot RMSF
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.plot(residue_ids, rmsf, 'b-', alpha=0.8, linewidth=1)
    ax.fill_between(residue_ids, 0, rmsf, alpha=0.3)
    ax.axhline(y=rmsf.mean(), color='r', linestyle='--', label=f'Mean: {rmsf.mean():.2f} Å')
    ax.axhline(y=threshold, color='orange', linestyle=':', label=f'Flexible threshold: {threshold:.2f} Å')
    
    ax.set_xlabel('Residue Number')
    ax.set_ylabel('RMSF (Å)')
    ax.set_title('Per-Residue RMSF - Model B (Arg-TRIS-PEG2) with GLUT1')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plot_path = os.path.join(OUTPUT_DIR, 'rmsf_analysis.png')
    plt.savefig(plot_path, dpi=150)
    print(f"\n✅ RMSF plot saved: {plot_path}")
    plt.close()
    
    # Also calculate ligand RMSF
    print("\n--- Ligand RMSF ---")
    ligand_atoms = traj.topology.select('not protein and not water and not (name Na or name Cl)')
    
    if len(ligand_atoms) > 0:
        # Superpose on protein, then calculate ligand RMSF
        traj.superpose(traj, frame=0, atom_indices=ca_atoms)
        
        lig_positions = traj.xyz[:, ligand_atoms, :]
        lig_mean = np.mean(lig_positions, axis=0)
        lig_fluct = lig_positions - lig_mean
        lig_rmsf = np.sqrt(np.mean(np.sum(lig_fluct**2, axis=2), axis=0)) * 10
        
        print(f"Ligand RMSF:")
        print(f"  Mean: {lig_rmsf.mean():.2f} Å")
        print(f"  Max: {lig_rmsf.max():.2f} Å")
        print(f"  Min: {lig_rmsf.min():.2f} Å")
    
    return rmsf, residue_ids


def main():
    print("=" * 60)
    print("H-bond and RMSF Analysis - Model B (Arg-TRIS-PEG2)")
    print("=" * 60)
    
    traj_path = os.path.join(RESULTS_DIR, "prod_model_b.dcd")
    topology_path = os.path.join(RESULTS_DIR, "prod_model_b_final.pdb")
    
    # 1. Hydrogen Bond Analysis
    hbond_counts, hbond_details = analyze_hydrogen_bonds(traj_path, topology_path)
    
    # 2. RMSF Analysis
    rmsf, residue_ids = analyze_rmsf(traj_path, topology_path)
    
    # Summary
    print("\n" + "=" * 60)
    print("Analysis Summary")
    print("=" * 60)
    
    print(f"\n📊 Hydrogen Bonds:")
    print(f"  Mean H-bonds per frame: {hbond_counts.mean():.1f} ± {hbond_counts.std():.1f}")
    print(f"  Stable H-bonds (>50% occupancy): {sum(1 for h in hbond_details if h['occupancy'] > 50)}")
    
    if rmsf is not None:
        print(f"\n📊 RMSF:")
        print(f"  Mean protein RMSF: {rmsf.mean():.2f} Å")
        print(f"  Most flexible residue: {rmsf.max():.2f} Å")
        print(f"  Most rigid residue: {rmsf.min():.2f} Å")
    
    # Save results
    results = {
        'model': 'Model B (Arg-TRIS-PEG2)',
        'hbond_mean': float(hbond_counts.mean()),
        'hbond_std': float(hbond_counts.std()),
        'hbond_details': hbond_details,
        'rmsf_mean': float(rmsf.mean()) if rmsf is not None else None,
        'rmsf_max': float(rmsf.max()) if rmsf is not None else None,
        'rmsf_min': float(rmsf.min()) if rmsf is not None else None
    }
    
    results_path = os.path.join(OUTPUT_DIR, 'hbond_rmsf_model_b.json')
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    print(f"\n✅ Results saved: {results_path}")


if __name__ == "__main__":
    main()
