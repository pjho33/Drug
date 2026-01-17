#!/usr/bin/env python
"""
Advanced Binding Analysis: Model A vs Model B
==============================================
Additional binding metrics:
1. Hydrogen Bonds (H-bonds)
2. RMSF (Root Mean Square Fluctuation) - ligand flexibility
3. Radius of Gyration - ligand compactness
4. Binding Site Residence Time
5. Salt Bridge / Electrostatic Interactions
6. Hydrophobic Contacts
7. Center of Mass Distance to Binding Site
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import json
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

# Paths
MODEL_A_DIR = "/home/pjho3/projects/Drug/results/phase2_gnina_md/model_a"
MODEL_B_DIR = "/home/pjho3/projects/Drug/results/phase2_gnina_md/model_b"
OUTPUT_DIR = "/home/pjho3/projects/Drug/results/analysis/comparison"

os.makedirs(OUTPUT_DIR, exist_ok=True)


def load_trajectory(results_dir, model_name):
    """Load trajectory with topology"""
    dcd_path = os.path.join(results_dir, f"prod_{model_name}.dcd")
    pdb_path = os.path.join(results_dir, f"prod_{model_name}_final.pdb")
    
    print(f"\nLoading {model_name}...")
    traj = md.load(dcd_path, top=pdb_path)
    print(f"  Frames: {traj.n_frames}, Atoms: {traj.n_atoms}")
    
    return traj


def analyze_hydrogen_bonds(traj, model_name):
    """Analyze protein-ligand hydrogen bonds"""
    print(f"\n{model_name} Hydrogen Bond Analysis:")
    
    # Get ligand and protein indices
    ligand_idx = traj.topology.select('not protein and not water and not (name Na or name Cl)')
    protein_idx = traj.topology.select('protein')
    
    # Find H-bond donors and acceptors in ligand
    ligand_donors = []
    ligand_acceptors = []
    
    for idx in ligand_idx:
        atom = traj.topology.atom(idx)
        # Donors: N-H, O-H
        if atom.element.symbol in ['N', 'O']:
            # Check if has hydrogen neighbor
            for neighbor in atom.residue.atoms:
                if neighbor.element.symbol == 'H':
                    ligand_donors.append(idx)
                    break
            ligand_acceptors.append(idx)
    
    print(f"  Ligand H-bond donors: {len(ligand_donors)}")
    print(f"  Ligand H-bond acceptors: {len(ligand_acceptors)}")
    
    # Use mdtraj's hydrogen bond detection
    # This finds all H-bonds in the system
    hbonds = md.baker_hubbard(traj, freq=0.1)  # At least 10% of frames
    
    # Filter for protein-ligand H-bonds
    ligand_set = set(ligand_idx)
    protein_set = set(protein_idx)
    
    pl_hbonds = []
    for hb in hbonds:
        donor, hydrogen, acceptor = hb
        donor_in_lig = donor in ligand_set
        acceptor_in_lig = acceptor in ligand_set
        donor_in_prot = donor in protein_set
        acceptor_in_prot = acceptor in protein_set
        
        # Protein-ligand H-bond
        if (donor_in_lig and acceptor_in_prot) or (donor_in_prot and acceptor_in_lig):
            pl_hbonds.append(hb)
    
    print(f"  Protein-Ligand H-bonds (>10% occupancy): {len(pl_hbonds)}")
    
    # Count H-bonds per frame
    hbonds_per_frame = []
    for frame_idx in range(0, traj.n_frames, 5):  # Sample every 5th frame
        frame_hbonds = md.baker_hubbard(traj[frame_idx], freq=0.0)
        
        count = 0
        for hb in frame_hbonds:
            donor, hydrogen, acceptor = hb
            donor_in_lig = donor in ligand_set
            acceptor_in_lig = acceptor in ligand_set
            donor_in_prot = donor in protein_set
            acceptor_in_prot = acceptor in protein_set
            
            if (donor_in_lig and acceptor_in_prot) or (donor_in_prot and acceptor_in_lig):
                count += 1
        
        hbonds_per_frame.append(count)
    
    hbonds_per_frame = np.array(hbonds_per_frame)
    print(f"  Mean H-bonds per frame: {np.mean(hbonds_per_frame):.2f} ± {np.std(hbonds_per_frame):.2f}")
    
    return {
        'persistent_hbonds': len(pl_hbonds),
        'hbonds_per_frame': hbonds_per_frame,
        'mean': np.mean(hbonds_per_frame),
        'std': np.std(hbonds_per_frame)
    }


def analyze_rmsf(traj, model_name):
    """Analyze ligand RMSF (flexibility)"""
    print(f"\n{model_name} RMSF Analysis (Ligand Flexibility):")
    
    ligand_idx = traj.topology.select('not protein and not water and not (name Na or name Cl)')
    protein_ca = traj.topology.select('protein and name CA')
    
    # Align trajectory to protein
    traj_aligned = traj.superpose(traj, frame=0, atom_indices=protein_ca)
    
    # Calculate RMSF for ligand atoms
    ligand_positions = traj_aligned.xyz[:, ligand_idx, :]
    mean_positions = np.mean(ligand_positions, axis=0)
    
    # RMSF per atom
    rmsf = np.sqrt(np.mean(np.sum((ligand_positions - mean_positions)**2, axis=2), axis=0)) * 10  # nm to Å
    
    # Overall ligand RMSF
    mean_rmsf = np.mean(rmsf)
    max_rmsf = np.max(rmsf)
    
    print(f"  Mean ligand RMSF: {mean_rmsf:.2f} Å")
    print(f"  Max ligand RMSF: {max_rmsf:.2f} Å")
    print(f"  (Lower RMSF = more rigid/stable binding)")
    
    return {
        'rmsf_per_atom': rmsf,
        'mean_rmsf': mean_rmsf,
        'max_rmsf': max_rmsf
    }


def analyze_radius_of_gyration(traj, model_name):
    """Analyze ligand radius of gyration (compactness)"""
    print(f"\n{model_name} Radius of Gyration Analysis:")
    
    ligand_idx = traj.topology.select('not protein and not water and not (name Na or name Cl)')
    
    # Calculate Rg for each frame
    rg_values = []
    
    for frame_idx in range(traj.n_frames):
        positions = traj.xyz[frame_idx, ligand_idx, :]
        com = np.mean(positions, axis=0)
        rg = np.sqrt(np.mean(np.sum((positions - com)**2, axis=1))) * 10  # nm to Å
        rg_values.append(rg)
    
    rg_values = np.array(rg_values)
    
    print(f"  Mean Rg: {np.mean(rg_values):.2f} ± {np.std(rg_values):.2f} Å")
    print(f"  (Stable Rg = ligand maintains conformation)")
    
    return {
        'rg_values': rg_values,
        'mean': np.mean(rg_values),
        'std': np.std(rg_values)
    }


def analyze_com_distance(traj, model_name):
    """Analyze distance from ligand COM to binding site center"""
    print(f"\n{model_name} Center of Mass Distance Analysis:")
    
    ligand_idx = traj.topology.select('not protein and not water and not (name Na or name Cl)')
    
    # Define binding site as protein atoms within 10Å of ligand in first frame
    first_frame_lig_pos = traj.xyz[0, ligand_idx, :]
    lig_com_first = np.mean(first_frame_lig_pos, axis=0)
    
    protein_idx = traj.topology.select('protein and name CA')
    first_frame_prot_pos = traj.xyz[0, protein_idx, :]
    
    # Find binding site CA atoms (within 1.5 nm of ligand COM)
    dists = np.linalg.norm(first_frame_prot_pos - lig_com_first, axis=1)
    binding_site_idx = protein_idx[dists < 1.5]  # 15 Å
    
    print(f"  Binding site CA atoms: {len(binding_site_idx)}")
    
    # Calculate COM distance over trajectory
    com_distances = []
    
    for frame_idx in range(traj.n_frames):
        lig_pos = traj.xyz[frame_idx, ligand_idx, :]
        lig_com = np.mean(lig_pos, axis=0)
        
        bs_pos = traj.xyz[frame_idx, binding_site_idx, :]
        bs_com = np.mean(bs_pos, axis=0)
        
        dist = np.linalg.norm(lig_com - bs_com) * 10  # nm to Å
        com_distances.append(dist)
    
    com_distances = np.array(com_distances)
    
    print(f"  Mean COM distance: {np.mean(com_distances):.2f} ± {np.std(com_distances):.2f} Å")
    print(f"  (Lower = ligand stays in binding site)")
    
    return {
        'com_distances': com_distances,
        'mean': np.mean(com_distances),
        'std': np.std(com_distances)
    }


def analyze_hydrophobic_contacts(traj, model_name, cutoff=0.4):
    """Analyze hydrophobic contacts"""
    print(f"\n{model_name} Hydrophobic Contact Analysis:")
    
    # Hydrophobic residues
    hydrophobic_res = ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO']
    
    ligand_idx = traj.topology.select('not protein and not water and not (name Na or name Cl)')
    
    # Find hydrophobic ligand atoms (carbons not bonded to N/O)
    ligand_carbons = []
    for idx in ligand_idx:
        atom = traj.topology.atom(idx)
        if atom.element.symbol == 'C':
            ligand_carbons.append(idx)
    
    # Find hydrophobic protein atoms
    hydrophobic_protein = []
    for residue in traj.topology.residues:
        if residue.name in hydrophobic_res:
            for atom in residue.atoms:
                if atom.element.symbol == 'C':
                    hydrophobic_protein.append(atom.index)
    
    print(f"  Ligand carbons: {len(ligand_carbons)}")
    print(f"  Hydrophobic protein carbons: {len(hydrophobic_protein)}")
    
    # Count hydrophobic contacts per frame
    contacts_per_frame = []
    
    for frame_idx in range(0, traj.n_frames, 5):
        frame_pos = traj.xyz[frame_idx]
        lig_pos = frame_pos[ligand_carbons]
        prot_pos = frame_pos[hydrophobic_protein]
        
        count = 0
        for lp in lig_pos:
            dists = np.linalg.norm(prot_pos - lp, axis=1)
            count += np.sum(dists < cutoff)
        
        contacts_per_frame.append(count)
    
    contacts_per_frame = np.array(contacts_per_frame)
    
    print(f"  Mean hydrophobic contacts: {np.mean(contacts_per_frame):.1f} ± {np.std(contacts_per_frame):.1f}")
    
    return {
        'contacts': contacts_per_frame,
        'mean': np.mean(contacts_per_frame),
        'std': np.std(contacts_per_frame)
    }


def analyze_salt_bridges(traj, model_name, cutoff=0.4):
    """Analyze salt bridges / electrostatic interactions"""
    print(f"\n{model_name} Salt Bridge / Electrostatic Analysis:")
    
    # Charged residues
    positive_res = ['ARG', 'LYS', 'HIS']
    negative_res = ['ASP', 'GLU']
    
    ligand_idx = traj.topology.select('not protein and not water and not (name Na or name Cl)')
    
    # Find charged atoms in ligand
    ligand_charged = []
    for idx in ligand_idx:
        atom = traj.topology.atom(idx)
        if atom.element.symbol in ['N', 'O']:  # Potential charged groups
            ligand_charged.append(idx)
    
    # Find charged protein atoms
    positive_atoms = []
    negative_atoms = []
    
    for residue in traj.topology.residues:
        if residue.name in positive_res:
            for atom in residue.atoms:
                if atom.name in ['NZ', 'NH1', 'NH2', 'NE', 'ND1', 'NE2']:  # Charged nitrogens
                    positive_atoms.append(atom.index)
        elif residue.name in negative_res:
            for atom in residue.atoms:
                if atom.name in ['OD1', 'OD2', 'OE1', 'OE2']:  # Carboxylate oxygens
                    negative_atoms.append(atom.index)
    
    print(f"  Ligand charged atoms: {len(ligand_charged)}")
    print(f"  Protein positive atoms: {len(positive_atoms)}")
    print(f"  Protein negative atoms: {len(negative_atoms)}")
    
    # Count salt bridges per frame
    salt_bridges_per_frame = []
    
    for frame_idx in range(0, traj.n_frames, 5):
        frame_pos = traj.xyz[frame_idx]
        lig_pos = frame_pos[ligand_charged]
        
        count = 0
        # Ligand to positive protein
        if len(positive_atoms) > 0:
            pos_prot = frame_pos[positive_atoms]
            for lp in lig_pos:
                dists = np.linalg.norm(pos_prot - lp, axis=1)
                count += np.sum(dists < cutoff)
        
        # Ligand to negative protein
        if len(negative_atoms) > 0:
            neg_prot = frame_pos[negative_atoms]
            for lp in lig_pos:
                dists = np.linalg.norm(neg_prot - lp, axis=1)
                count += np.sum(dists < cutoff)
        
        salt_bridges_per_frame.append(count)
    
    salt_bridges_per_frame = np.array(salt_bridges_per_frame)
    
    print(f"  Mean electrostatic contacts: {np.mean(salt_bridges_per_frame):.1f} ± {np.std(salt_bridges_per_frame):.1f}")
    
    return {
        'contacts': salt_bridges_per_frame,
        'mean': np.mean(salt_bridges_per_frame),
        'std': np.std(salt_bridges_per_frame)
    }


def analyze_binding_stability(traj, model_name):
    """Analyze binding stability via ligand displacement over time"""
    print(f"\n{model_name} Binding Stability Analysis:")
    
    ligand_idx = traj.topology.select('not protein and not water and not (name Na or name Cl)')
    protein_ca = traj.topology.select('protein and name CA')
    
    # Align to protein
    traj_aligned = traj.superpose(traj, frame=0, atom_indices=protein_ca)
    
    # Calculate ligand COM displacement from initial position
    initial_com = np.mean(traj_aligned.xyz[0, ligand_idx, :], axis=0)
    
    displacements = []
    for frame_idx in range(traj.n_frames):
        current_com = np.mean(traj_aligned.xyz[frame_idx, ligand_idx, :], axis=0)
        disp = np.linalg.norm(current_com - initial_com) * 10  # nm to Å
        displacements.append(disp)
    
    displacements = np.array(displacements)
    
    # Calculate "escape" events (displacement > 5Å)
    escape_threshold = 5.0
    escape_frames = np.sum(displacements > escape_threshold)
    escape_fraction = escape_frames / len(displacements) * 100
    
    print(f"  Mean displacement: {np.mean(displacements):.2f} ± {np.std(displacements):.2f} Å")
    print(f"  Max displacement: {np.max(displacements):.2f} Å")
    print(f"  Escape events (>5Å): {escape_frames}/{len(displacements)} ({escape_fraction:.1f}%)")
    
    return {
        'displacements': displacements,
        'mean': np.mean(displacements),
        'std': np.std(displacements),
        'max': np.max(displacements),
        'escape_fraction': escape_fraction
    }


def plot_advanced_comparison(results_a, results_b):
    """Create advanced comparison plots"""
    fig, axes = plt.subplots(3, 3, figsize=(16, 14))
    
    # 1. H-bonds over time
    ax = axes[0, 0]
    time_a = np.arange(len(results_a['hbonds']['hbonds_per_frame'])) * 0.05
    time_b = np.arange(len(results_b['hbonds']['hbonds_per_frame'])) * 0.05
    ax.plot(time_a, results_a['hbonds']['hbonds_per_frame'], 'b-', alpha=0.7, label='Model A')
    ax.plot(time_b, results_b['hbonds']['hbonds_per_frame'], 'r-', alpha=0.7, label='Model B')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('H-bonds')
    ax.set_title('Protein-Ligand Hydrogen Bonds')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 2. RMSF comparison (bar chart)
    ax = axes[0, 1]
    models = ['Model A\n(TRIS)', 'Model B\n(Arg-TRIS)']
    rmsf_vals = [results_a['rmsf']['mean_rmsf'], results_b['rmsf']['mean_rmsf']]
    colors = ['blue', 'red']
    bars = ax.bar(models, rmsf_vals, color=colors, alpha=0.7)
    ax.set_ylabel('Mean RMSF (Å)')
    ax.set_title('Ligand Flexibility (RMSF)')
    for bar, val in zip(bars, rmsf_vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05, f'{val:.2f}', ha='center')
    
    # 3. Radius of Gyration over time
    ax = axes[0, 2]
    time_a = np.arange(len(results_a['rg']['rg_values'])) * 0.01
    time_b = np.arange(len(results_b['rg']['rg_values'])) * 0.01
    ax.plot(time_a, results_a['rg']['rg_values'], 'b-', alpha=0.7, label='Model A')
    ax.plot(time_b, results_b['rg']['rg_values'], 'r-', alpha=0.7, label='Model B')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Radius of Gyration (Å)')
    ax.set_title('Ligand Compactness')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 4. COM distance over time
    ax = axes[1, 0]
    time_a = np.arange(len(results_a['com']['com_distances'])) * 0.01
    time_b = np.arange(len(results_b['com']['com_distances'])) * 0.01
    ax.plot(time_a, results_a['com']['com_distances'], 'b-', alpha=0.7, label='Model A')
    ax.plot(time_b, results_b['com']['com_distances'], 'r-', alpha=0.7, label='Model B')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('COM Distance (Å)')
    ax.set_title('Distance to Binding Site Center')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 5. Hydrophobic contacts
    ax = axes[1, 1]
    time_a = np.arange(len(results_a['hydrophobic']['contacts'])) * 0.05
    time_b = np.arange(len(results_b['hydrophobic']['contacts'])) * 0.05
    ax.plot(time_a, results_a['hydrophobic']['contacts'], 'b-', alpha=0.7, label='Model A')
    ax.plot(time_b, results_b['hydrophobic']['contacts'], 'r-', alpha=0.7, label='Model B')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Hydrophobic Contacts')
    ax.set_title('Hydrophobic Interactions')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 6. Salt bridges / electrostatic
    ax = axes[1, 2]
    time_a = np.arange(len(results_a['salt_bridges']['contacts'])) * 0.05
    time_b = np.arange(len(results_b['salt_bridges']['contacts'])) * 0.05
    ax.plot(time_a, results_a['salt_bridges']['contacts'], 'b-', alpha=0.7, label='Model A')
    ax.plot(time_b, results_b['salt_bridges']['contacts'], 'r-', alpha=0.7, label='Model B')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Electrostatic Contacts')
    ax.set_title('Salt Bridges / Electrostatic')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 7. Binding stability (displacement)
    ax = axes[2, 0]
    time_a = np.arange(len(results_a['stability']['displacements'])) * 0.01
    time_b = np.arange(len(results_b['stability']['displacements'])) * 0.01
    ax.plot(time_a, results_a['stability']['displacements'], 'b-', alpha=0.7, label='Model A')
    ax.plot(time_b, results_b['stability']['displacements'], 'r-', alpha=0.7, label='Model B')
    ax.axhline(y=5.0, color='k', linestyle='--', alpha=0.5, label='Escape threshold')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Displacement (Å)')
    ax.set_title('Ligand Displacement from Initial Position')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 8. Summary comparison (spider/radar chart alternative - bar chart)
    ax = axes[2, 1]
    metrics = ['H-bonds', 'RMSF\n(inv)', 'Rg\nstability', 'Hydrophobic', 'Electrostatic', 'Stability\n(inv)']
    
    # Normalize values (higher is better)
    model_a_vals = [
        results_a['hbonds']['mean'],
        1.0 / (results_a['rmsf']['mean_rmsf'] + 0.1),  # Inverse - lower is better
        1.0 / (results_a['rg']['std'] + 0.1),  # Inverse of std - lower is better
        results_a['hydrophobic']['mean'],
        results_a['salt_bridges']['mean'],
        1.0 / (results_a['stability']['mean'] + 0.1)  # Inverse - lower is better
    ]
    model_b_vals = [
        results_b['hbonds']['mean'],
        1.0 / (results_b['rmsf']['mean_rmsf'] + 0.1),
        1.0 / (results_b['rg']['std'] + 0.1),
        results_b['hydrophobic']['mean'],
        results_b['salt_bridges']['mean'],
        1.0 / (results_b['stability']['mean'] + 0.1)
    ]
    
    # Normalize to 0-100 scale
    max_vals = [max(a, b) for a, b in zip(model_a_vals, model_b_vals)]
    model_a_norm = [v/m*100 if m > 0 else 0 for v, m in zip(model_a_vals, max_vals)]
    model_b_norm = [v/m*100 if m > 0 else 0 for v, m in zip(model_b_vals, max_vals)]
    
    x = np.arange(len(metrics))
    width = 0.35
    ax.bar(x - width/2, model_a_norm, width, label='Model A', color='blue', alpha=0.7)
    ax.bar(x + width/2, model_b_norm, width, label='Model B', color='red', alpha=0.7)
    ax.set_ylabel('Normalized Score')
    ax.set_title('Binding Quality Metrics (higher = better)')
    ax.set_xticks(x)
    ax.set_xticklabels(metrics, fontsize=8)
    ax.legend()
    ax.set_ylim(0, 120)
    
    # 9. Overall binding score
    ax = axes[2, 2]
    
    # Calculate composite binding score
    def calc_binding_score(r):
        score = 0
        score += r['hbonds']['mean'] * 10  # H-bonds are important
        score += 100 / (r['rmsf']['mean_rmsf'] + 1)  # Lower RMSF is better
        score += r['hydrophobic']['mean'] * 0.5
        score += r['salt_bridges']['mean'] * 5
        score += 100 / (r['stability']['mean'] + 1)  # Lower displacement is better
        return score
    
    score_a = calc_binding_score(results_a)
    score_b = calc_binding_score(results_b)
    
    bars = ax.bar(['Model A\n(TRIS)', 'Model B\n(Arg-TRIS)'], [score_a, score_b], 
                  color=['blue', 'red'], alpha=0.7)
    ax.set_ylabel('Composite Binding Score')
    ax.set_title('Overall Binding Quality')
    for bar, val in zip(bars, [score_a, score_b]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5, f'{val:.1f}', ha='center', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'advanced_binding_comparison.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n📊 Advanced comparison plot saved!")
    
    return score_a, score_b


def main():
    print("=" * 70)
    print("Advanced Binding Analysis: Model A vs Model B")
    print("=" * 70)
    
    # Load trajectories
    traj_a = load_trajectory(MODEL_A_DIR, "model_a")
    traj_b = load_trajectory(MODEL_B_DIR, "model_b")
    
    # Analyze Model A
    print("\n" + "=" * 50)
    print("ANALYZING MODEL A (Basic TRIS)")
    print("=" * 50)
    
    results_a = {
        'hbonds': analyze_hydrogen_bonds(traj_a, "Model A"),
        'rmsf': analyze_rmsf(traj_a, "Model A"),
        'rg': analyze_radius_of_gyration(traj_a, "Model A"),
        'com': analyze_com_distance(traj_a, "Model A"),
        'hydrophobic': analyze_hydrophobic_contacts(traj_a, "Model A"),
        'salt_bridges': analyze_salt_bridges(traj_a, "Model A"),
        'stability': analyze_binding_stability(traj_a, "Model A")
    }
    
    # Analyze Model B
    print("\n" + "=" * 50)
    print("ANALYZING MODEL B (Arg-TRIS)")
    print("=" * 50)
    
    results_b = {
        'hbonds': analyze_hydrogen_bonds(traj_b, "Model B"),
        'rmsf': analyze_rmsf(traj_b, "Model B"),
        'rg': analyze_radius_of_gyration(traj_b, "Model B"),
        'com': analyze_com_distance(traj_b, "Model B"),
        'hydrophobic': analyze_hydrophobic_contacts(traj_b, "Model B"),
        'salt_bridges': analyze_salt_bridges(traj_b, "Model B"),
        'stability': analyze_binding_stability(traj_b, "Model B")
    }
    
    # Plot comparison
    score_a, score_b = plot_advanced_comparison(results_a, results_b)
    
    # Summary table
    print("\n" + "=" * 80)
    print("ADVANCED BINDING METRICS COMPARISON")
    print("=" * 80)
    
    print("\n┌────────────────────────────┬───────────────────┬───────────────────┬─────────┐")
    print("│ Metric                     │ Model A (TRIS)    │ Model B (Arg)     │ Winner  │")
    print("├────────────────────────────┼───────────────────┼───────────────────┼─────────┤")
    
    # H-bonds
    winner = "B ✓" if results_b['hbonds']['mean'] > results_a['hbonds']['mean'] else "A ✓"
    print(f"│ H-bonds (mean)             │ {results_a['hbonds']['mean']:6.2f} ± {results_a['hbonds']['std']:5.2f}   │ {results_b['hbonds']['mean']:6.2f} ± {results_b['hbonds']['std']:5.2f}   │   {winner}   │")
    
    # RMSF (lower is better)
    winner = "B ✓" if results_b['rmsf']['mean_rmsf'] < results_a['rmsf']['mean_rmsf'] else "A ✓"
    print(f"│ RMSF (Å, lower=stable)     │ {results_a['rmsf']['mean_rmsf']:6.2f}            │ {results_b['rmsf']['mean_rmsf']:6.2f}            │   {winner}   │")
    
    # Rg stability (lower std is better)
    winner = "B ✓" if results_b['rg']['std'] < results_a['rg']['std'] else "A ✓"
    print(f"│ Rg std (Å, lower=stable)   │ {results_a['rg']['std']:6.2f}            │ {results_b['rg']['std']:6.2f}            │   {winner}   │")
    
    # COM distance (lower is better)
    winner = "B ✓" if results_b['com']['mean'] < results_a['com']['mean'] else "A ✓"
    print(f"│ COM dist (Å, lower=bound)  │ {results_a['com']['mean']:6.2f} ± {results_a['com']['std']:5.2f}   │ {results_b['com']['mean']:6.2f} ± {results_b['com']['std']:5.2f}   │   {winner}   │")
    
    # Hydrophobic contacts
    winner = "B ✓" if results_b['hydrophobic']['mean'] > results_a['hydrophobic']['mean'] else "A ✓"
    print(f"│ Hydrophobic contacts       │ {results_a['hydrophobic']['mean']:6.1f} ± {results_a['hydrophobic']['std']:5.1f}   │ {results_b['hydrophobic']['mean']:6.1f} ± {results_b['hydrophobic']['std']:5.1f}   │   {winner}   │")
    
    # Salt bridges
    winner = "B ✓" if results_b['salt_bridges']['mean'] > results_a['salt_bridges']['mean'] else "A ✓"
    print(f"│ Electrostatic contacts     │ {results_a['salt_bridges']['mean']:6.1f} ± {results_a['salt_bridges']['std']:5.1f}   │ {results_b['salt_bridges']['mean']:6.1f} ± {results_b['salt_bridges']['std']:5.1f}   │   {winner}   │")
    
    # Stability (lower displacement is better)
    winner = "B ✓" if results_b['stability']['mean'] < results_a['stability']['mean'] else "A ✓"
    print(f"│ Displacement (Å, lower=ok) │ {results_a['stability']['mean']:6.2f} ± {results_a['stability']['std']:5.2f}   │ {results_b['stability']['mean']:6.2f} ± {results_b['stability']['std']:5.2f}   │   {winner}   │")
    
    print("├────────────────────────────┼───────────────────┼───────────────────┼─────────┤")
    winner = "B ✓" if score_b > score_a else "A ✓"
    print(f"│ COMPOSITE SCORE            │ {score_a:6.1f}            │ {score_b:6.1f}            │   {winner}   │")
    print("└────────────────────────────┴───────────────────┴───────────────────┴─────────┘")
    
    # Save results
    summary = {
        'model_a': {
            'hbonds_mean': float(results_a['hbonds']['mean']),
            'rmsf_mean': float(results_a['rmsf']['mean_rmsf']),
            'rg_mean': float(results_a['rg']['mean']),
            'rg_std': float(results_a['rg']['std']),
            'com_mean': float(results_a['com']['mean']),
            'hydrophobic_mean': float(results_a['hydrophobic']['mean']),
            'salt_bridges_mean': float(results_a['salt_bridges']['mean']),
            'displacement_mean': float(results_a['stability']['mean']),
            'composite_score': float(score_a)
        },
        'model_b': {
            'hbonds_mean': float(results_b['hbonds']['mean']),
            'rmsf_mean': float(results_b['rmsf']['mean_rmsf']),
            'rg_mean': float(results_b['rg']['mean']),
            'rg_std': float(results_b['rg']['std']),
            'com_mean': float(results_b['com']['mean']),
            'hydrophobic_mean': float(results_b['hydrophobic']['mean']),
            'salt_bridges_mean': float(results_b['salt_bridges']['mean']),
            'displacement_mean': float(results_b['stability']['mean']),
            'composite_score': float(score_b)
        }
    }
    
    with open(os.path.join(OUTPUT_DIR, 'advanced_binding_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n💾 Results saved to: {OUTPUT_DIR}/advanced_binding_summary.json")
    print(f"📊 Plot saved to: {OUTPUT_DIR}/advanced_binding_comparison.png")
    
    print("\n" + "=" * 70)
    print("Analysis Complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
