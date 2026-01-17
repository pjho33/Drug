#!/usr/bin/env python
"""
Model A vs Model B Comparison Analysis
=======================================
Compares TRIS (Model A) vs Arg-TRIS (Model B) MD simulations:
- RMSD analysis
- Binding energy estimation
- Cation-π interactions (Model B only)
- Contact analysis
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import json
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
    # Use final PDB which has the solvated system topology
    pdb_path = os.path.join(results_dir, f"prod_{model_name}_final.pdb")
    
    print(f"\nLoading {model_name}...")
    print(f"  DCD: {dcd_path}")
    print(f"  PDB: {pdb_path}")
    
    traj = md.load(dcd_path, top=pdb_path)
    print(f"  Frames: {traj.n_frames}")
    print(f"  Atoms: {traj.n_atoms}")
    
    return traj


def analyze_rmsd(traj, model_name):
    """Calculate RMSD for ligand and protein"""
    # Select ligand atoms
    ligand_idx = traj.topology.select('not protein and not water and not (name Na or name Cl)')
    protein_idx = traj.topology.select('protein and name CA')
    
    print(f"\n{model_name} RMSD Analysis:")
    print(f"  Ligand atoms: {len(ligand_idx)}")
    print(f"  Protein CA atoms: {len(protein_idx)}")
    
    # Align to first frame using protein CA
    traj_aligned = traj.superpose(traj, frame=0, atom_indices=protein_idx)
    
    # Calculate ligand RMSD
    ligand_rmsd = md.rmsd(traj_aligned, traj_aligned, frame=0, atom_indices=ligand_idx) * 10  # nm to Å
    
    # Calculate protein RMSD
    protein_rmsd = md.rmsd(traj_aligned, traj_aligned, frame=0, atom_indices=protein_idx) * 10  # nm to Å
    
    print(f"  Ligand RMSD: {np.mean(ligand_rmsd):.2f} ± {np.std(ligand_rmsd):.2f} Å")
    print(f"  Protein RMSD: {np.mean(protein_rmsd):.2f} ± {np.std(protein_rmsd):.2f} Å")
    
    return {
        'ligand_rmsd': ligand_rmsd,
        'protein_rmsd': protein_rmsd,
        'ligand_mean': np.mean(ligand_rmsd),
        'ligand_std': np.std(ligand_rmsd),
        'protein_mean': np.mean(protein_rmsd),
        'protein_std': np.std(protein_rmsd)
    }


def analyze_contacts(traj, model_name, cutoff=0.5):
    """Analyze protein-ligand contacts using mdtraj's efficient compute_neighbors"""
    ligand_idx = traj.topology.select('not protein and not water and not (name Na or name Cl)')
    protein_idx = traj.topology.select('protein and element != H')
    
    print(f"\n{model_name} Contact Analysis (cutoff={cutoff*10:.1f} Å):")
    print(f"  Ligand atoms: {len(ligand_idx)}, Protein heavy atoms: {len(protein_idx)}")
    
    # Use mdtraj's efficient contact calculation
    contacts_per_frame = []
    
    # Sample frames for efficiency (every 5th frame)
    frame_indices = list(range(0, traj.n_frames, 5))
    
    for i, frame_idx in enumerate(frame_indices):
        frame_pos = traj.xyz[frame_idx]
        ligand_pos = frame_pos[ligand_idx]
        protein_pos = frame_pos[protein_idx]
        
        # Vectorized distance calculation
        n_contacts = 0
        for lig_pos in ligand_pos:
            dists = np.linalg.norm(protein_pos - lig_pos, axis=1)
            n_contacts += np.sum(dists < cutoff)
        
        contacts_per_frame.append(n_contacts)
        
        if (i + 1) % 5 == 0:
            print(f"    Contact analysis: {(i+1)/len(frame_indices)*100:.0f}%")
    
    contacts = np.array(contacts_per_frame)
    
    print(f"  Mean contacts: {np.mean(contacts):.1f} ± {np.std(contacts):.1f}")
    
    return {
        'contacts': contacts,
        'mean': np.mean(contacts),
        'std': np.std(contacts)
    }


def analyze_cation_pi(traj, model_name):
    """Analyze cation-π interactions (for Model B with Arginine)"""
    from collections import Counter
    
    # Find nitrogen atoms in ligand (potential guanidinium)
    ligand_idx = traj.topology.select('not protein and not water and not (name Na or name Cl)')
    
    ligand_nitrogens = []
    for idx in ligand_idx:
        atom = traj.topology.atom(idx)
        if atom.element.symbol == 'N':
            ligand_nitrogens.append(idx)
    
    if len(ligand_nitrogens) == 0:
        print(f"\n{model_name}: No nitrogen atoms found in ligand (no cation-π possible)")
        return None
    
    # Find aromatic residues - limit to first 100 residues near binding site
    aromatic_residues = []
    for residue in traj.topology.residues:
        if residue.name in ['TRP', 'PHE', 'TYR']:
            ring_atoms = []
            for atom in residue.atoms:
                if residue.name == 'TRP':
                    if atom.name in ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']:
                        ring_atoms.append(atom.index)
                else:
                    if atom.name in ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']:
                        ring_atoms.append(atom.index)
            
            if len(ring_atoms) > 0:
                aromatic_residues.append({
                    'residue': f"{residue.name}{residue.resSeq}",
                    'ring_atoms': ring_atoms
                })
    
    print(f"\n{model_name} Cation-π Analysis:")
    print(f"  Ligand nitrogens: {len(ligand_nitrogens)}")
    print(f"  Aromatic residues: {len(aromatic_residues)}")
    
    # Calculate minimum distance to aromatic rings per frame (sample every 5th)
    min_distances = []
    closest_residues = []
    
    frame_indices = list(range(0, traj.n_frames, 5))
    
    for frame_idx in frame_indices:
        frame_pos = traj.xyz[frame_idx]
        
        min_dist = float('inf')
        closest_res = None
        
        for n_idx in ligand_nitrogens:
            n_pos = frame_pos[n_idx]
            
            for arom in aromatic_residues:
                ring_center = np.mean(frame_pos[arom['ring_atoms']], axis=0)
                dist = np.linalg.norm(n_pos - ring_center)
                
                if dist < min_dist:
                    min_dist = dist
                    closest_res = arom['residue']
        
        min_distances.append(min_dist * 10)  # nm to Å
        closest_residues.append(closest_res)
    
    min_distances = np.array(min_distances)
    
    # Count frames with cation-π interaction (< 5 Å)
    cation_pi_frames = np.sum(min_distances < 5.0)
    cation_pi_fraction = cation_pi_frames / len(min_distances) * 100
    
    print(f"  Mean distance to aromatic: {np.mean(min_distances):.2f} ± {np.std(min_distances):.2f} Å")
    print(f"  Cation-π frames (<5Å): {cation_pi_frames}/{len(min_distances)} ({cation_pi_fraction:.1f}%)")
    
    # Most common closest residue
    residue_counts = Counter(closest_residues)
    print(f"  Most common closest aromatic: {residue_counts.most_common(3)}")
    
    return {
        'min_distances': min_distances,
        'mean_distance': np.mean(min_distances),
        'std_distance': np.std(min_distances),
        'cation_pi_fraction': cation_pi_fraction,
        'closest_residues': residue_counts.most_common(5)
    }


def estimate_binding_energy(traj, model_name):
    """
    Estimate relative binding energy using contact-based scoring
    This is a simplified proxy - not actual MM/PBSA
    """
    ligand_idx = traj.topology.select('not protein and not water and not (name Na or name Cl)')
    
    # Get heavy atoms only
    ligand_heavy = [i for i in ligand_idx if traj.topology.atom(i).element.symbol != 'H']
    
    # Find protein atoms within 6Å of ligand in any frame
    protein_idx = traj.topology.select('protein')
    protein_heavy = [i for i in protein_idx if traj.topology.atom(i).element.symbol != 'H']
    
    # Calculate average interaction score
    scores = []
    
    for frame_idx in range(0, traj.n_frames, 10):  # Sample every 10th frame
        frame_pos = traj.xyz[frame_idx]
        ligand_pos = frame_pos[ligand_heavy]
        protein_pos = frame_pos[protein_heavy]
        
        # Simple distance-based scoring
        score = 0
        for lig_pos in ligand_pos:
            dists = np.linalg.norm(protein_pos - lig_pos, axis=1)
            # Count favorable contacts (3-5 Å range)
            favorable = np.sum((dists > 0.3) & (dists < 0.5))
            # Penalize clashes (< 2.5 Å)
            clashes = np.sum(dists < 0.25)
            score += favorable - 2 * clashes
        
        scores.append(score)
    
    scores = np.array(scores)
    
    print(f"\n{model_name} Binding Score (contact-based proxy):")
    print(f"  Mean score: {np.mean(scores):.1f} ± {np.std(scores):.1f}")
    
    return {
        'scores': scores,
        'mean': np.mean(scores),
        'std': np.std(scores)
    }


def plot_comparison(results_a, results_b):
    """Create comparison plots"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    time_a = np.arange(len(results_a['rmsd']['ligand_rmsd'])) * 0.01  # Assuming 10 ps save interval
    time_b = np.arange(len(results_b['rmsd']['ligand_rmsd'])) * 0.01
    
    # 1. Ligand RMSD comparison
    ax1 = axes[0, 0]
    ax1.plot(time_a, results_a['rmsd']['ligand_rmsd'], 'b-', alpha=0.7, label='Model A (Basic TRIS)')
    ax1.plot(time_b, results_b['rmsd']['ligand_rmsd'], 'r-', alpha=0.7, label='Model B (Arg-TRIS)')
    ax1.set_xlabel('Time (ns)')
    ax1.set_ylabel('Ligand RMSD (Å)')
    ax1.set_title('Ligand RMSD Comparison')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Contact analysis (use sampled time)
    ax2 = axes[0, 1]
    time_contacts_a = np.arange(len(results_a['contacts']['contacts'])) * 0.05  # Sampled every 5 frames
    time_contacts_b = np.arange(len(results_b['contacts']['contacts'])) * 0.05
    ax2.plot(time_contacts_a, results_a['contacts']['contacts'], 'b-', alpha=0.7, label='Model A')
    ax2.plot(time_contacts_b, results_b['contacts']['contacts'], 'r-', alpha=0.7, label='Model B')
    ax2.set_xlabel('Time (ns)')
    ax2.set_ylabel('Number of Contacts')
    ax2.set_title('Protein-Ligand Contacts')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Cation-π distance (Model B only, or comparison if both have)
    ax3 = axes[1, 0]
    if results_b['cation_pi'] is not None:
        time_cation_b = np.arange(len(results_b['cation_pi']['min_distances'])) * 0.05  # Sampled every 5 frames
        ax3.plot(time_cation_b, results_b['cation_pi']['min_distances'], 'r-', alpha=0.7, label='Model B (Arg-TRIS)')
        ax3.axhline(y=5.0, color='k', linestyle='--', alpha=0.5, label='Cation-π cutoff (5Å)')
        ax3.set_xlabel('Time (ns)')
        ax3.set_ylabel('Distance to Aromatic (Å)')
        ax3.set_title('Cation-π Interaction Distance (Model B)')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
    else:
        ax3.text(0.5, 0.5, 'No cation-π data', ha='center', va='center', transform=ax3.transAxes)
    
    # 4. Summary bar chart
    ax4 = axes[1, 1]
    metrics = ['Ligand RMSD\n(Å)', 'Contacts', 'Binding Score']
    model_a_vals = [
        results_a['rmsd']['ligand_mean'],
        results_a['contacts']['mean'],
        results_a['binding']['mean']
    ]
    model_b_vals = [
        results_b['rmsd']['ligand_mean'],
        results_b['contacts']['mean'],
        results_b['binding']['mean']
    ]
    
    x = np.arange(len(metrics))
    width = 0.35
    
    # Normalize for visualization
    max_vals = [max(a, b) for a, b in zip(model_a_vals, model_b_vals)]
    model_a_norm = [v/m*100 if m > 0 else 0 for v, m in zip(model_a_vals, max_vals)]
    model_b_norm = [v/m*100 if m > 0 else 0 for v, m in zip(model_b_vals, max_vals)]
    
    bars1 = ax4.bar(x - width/2, model_a_norm, width, label='Model A (Basic TRIS)', color='blue', alpha=0.7)
    bars2 = ax4.bar(x + width/2, model_b_norm, width, label='Model B (Arg-TRIS)', color='red', alpha=0.7)
    
    # Add value labels
    for bar, val in zip(bars1, model_a_vals):
        ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2, f'{val:.1f}', 
                ha='center', va='bottom', fontsize=9)
    for bar, val in zip(bars2, model_b_vals):
        ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2, f'{val:.1f}', 
                ha='center', va='bottom', fontsize=9)
    
    ax4.set_ylabel('Normalized Value (%)')
    ax4.set_title('Summary Comparison')
    ax4.set_xticks(x)
    ax4.set_xticklabels(metrics)
    ax4.legend()
    ax4.set_ylim(0, 130)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'model_a_vs_b_comparison.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n📊 Comparison plot saved to: {OUTPUT_DIR}/model_a_vs_b_comparison.png")


def main():
    print("=" * 70)
    print("Model A vs Model B Comparison Analysis")
    print("Model A: Basic TRIS-PEG2 (no Arginine)")
    print("Model B: Arg-TRIS-PEG2 (with Arginine - Cation-π)")
    print("=" * 70)
    
    # Load trajectories
    traj_a = load_trajectory(MODEL_A_DIR, "model_a")
    traj_b = load_trajectory(MODEL_B_DIR, "model_b")
    
    # Analyze Model A
    print("\n" + "=" * 50)
    print("ANALYZING MODEL A (Basic TRIS)")
    print("=" * 50)
    
    rmsd_a = analyze_rmsd(traj_a, "Model A")
    contacts_a = analyze_contacts(traj_a, "Model A")
    cation_pi_a = analyze_cation_pi(traj_a, "Model A")
    binding_a = estimate_binding_energy(traj_a, "Model A")
    
    results_a = {
        'rmsd': rmsd_a,
        'contacts': contacts_a,
        'cation_pi': cation_pi_a,
        'binding': binding_a
    }
    
    # Analyze Model B
    print("\n" + "=" * 50)
    print("ANALYZING MODEL B (Arg-TRIS)")
    print("=" * 50)
    
    rmsd_b = analyze_rmsd(traj_b, "Model B")
    contacts_b = analyze_contacts(traj_b, "Model B")
    cation_pi_b = analyze_cation_pi(traj_b, "Model B")
    binding_b = estimate_binding_energy(traj_b, "Model B")
    
    results_b = {
        'rmsd': rmsd_b,
        'contacts': contacts_b,
        'cation_pi': cation_pi_b,
        'binding': binding_b
    }
    
    # Create comparison plots
    plot_comparison(results_a, results_b)
    
    # Summary report
    print("\n" + "=" * 70)
    print("COMPARISON SUMMARY")
    print("=" * 70)
    
    print("\n┌─────────────────────────────────────────────────────────────────┐")
    print("│                    Model A vs Model B Summary                    │")
    print("├─────────────────────┬──────────────────┬──────────────────┬──────┤")
    print("│ Metric              │ Model A (TRIS)   │ Model B (Arg)    │ Diff │")
    print("├─────────────────────┼──────────────────┼──────────────────┼──────┤")
    
    # RMSD
    rmsd_diff = rmsd_b['ligand_mean'] - rmsd_a['ligand_mean']
    print(f"│ Ligand RMSD (Å)     │ {rmsd_a['ligand_mean']:6.2f} ± {rmsd_a['ligand_std']:5.2f}  │ {rmsd_b['ligand_mean']:6.2f} ± {rmsd_b['ligand_std']:5.2f}  │{rmsd_diff:+5.1f} │")
    
    # Contacts
    contact_diff = contacts_b['mean'] - contacts_a['mean']
    print(f"│ Contacts            │ {contacts_a['mean']:6.1f} ± {contacts_a['std']:5.1f}  │ {contacts_b['mean']:6.1f} ± {contacts_b['std']:5.1f}  │{contact_diff:+5.0f} │")
    
    # Binding score
    binding_diff = binding_b['mean'] - binding_a['mean']
    print(f"│ Binding Score       │ {binding_a['mean']:6.1f} ± {binding_a['std']:5.1f}  │ {binding_b['mean']:6.1f} ± {binding_b['std']:5.1f}  │{binding_diff:+5.0f} │")
    
    # Cation-π (Model B only)
    if cation_pi_b is not None:
        print(f"│ Cation-π (<5Å)      │       N/A        │ {cation_pi_b['cation_pi_fraction']:6.1f}%          │  -   │")
    
    print("└─────────────────────┴──────────────────┴──────────────────┴──────┘")
    
    # Interpretation
    print("\n📋 INTERPRETATION:")
    print("-" * 50)
    
    if rmsd_b['ligand_mean'] < rmsd_a['ligand_mean']:
        print("✅ Model B shows LOWER ligand RMSD → more stable binding")
    else:
        print("⚠️  Model A shows lower ligand RMSD")
    
    if contacts_b['mean'] > contacts_a['mean']:
        print("✅ Model B has MORE protein-ligand contacts")
    else:
        print("⚠️  Model A has more contacts")
    
    if binding_b['mean'] > binding_a['mean']:
        print("✅ Model B has HIGHER binding score")
    else:
        print("⚠️  Model A has higher binding score")
    
    if cation_pi_b is not None and cation_pi_b['cation_pi_fraction'] > 20:
        print(f"✅ Model B shows significant cation-π interactions ({cation_pi_b['cation_pi_fraction']:.1f}% of frames)")
    
    # Save results
    summary = {
        'model_a': {
            'ligand_rmsd_mean': float(rmsd_a['ligand_mean']),
            'ligand_rmsd_std': float(rmsd_a['ligand_std']),
            'contacts_mean': float(contacts_a['mean']),
            'contacts_std': float(contacts_a['std']),
            'binding_score_mean': float(binding_a['mean']),
            'binding_score_std': float(binding_a['std'])
        },
        'model_b': {
            'ligand_rmsd_mean': float(rmsd_b['ligand_mean']),
            'ligand_rmsd_std': float(rmsd_b['ligand_std']),
            'contacts_mean': float(contacts_b['mean']),
            'contacts_std': float(contacts_b['std']),
            'binding_score_mean': float(binding_b['mean']),
            'binding_score_std': float(binding_b['std']),
            'cation_pi_fraction': float(cation_pi_b['cation_pi_fraction']) if cation_pi_b else None,
            'cation_pi_distance_mean': float(cation_pi_b['mean_distance']) if cation_pi_b else None
        },
        'comparison': {
            'rmsd_diff': float(rmsd_diff),
            'contacts_diff': float(contact_diff),
            'binding_diff': float(binding_diff)
        }
    }
    
    with open(os.path.join(OUTPUT_DIR, 'comparison_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n💾 Results saved to: {OUTPUT_DIR}/comparison_summary.json")
    print("\n" + "=" * 70)
    print("Analysis Complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
