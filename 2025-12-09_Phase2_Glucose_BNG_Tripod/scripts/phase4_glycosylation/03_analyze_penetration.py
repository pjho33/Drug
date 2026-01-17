#!/usr/bin/env python
"""
Phase 4: Analyze Glycan Penetration Results
============================================
Compare drug behavior between naked and glycosylated GLUT1:
1. COM distance to binding site over time
2. H-bonds with glycan vs protein
3. Solvent Accessible Surface Area (SAS)
4. Contact analysis
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
OUTPUT_DIR = "/home/pjho3/projects/Drug/results/phase4_glycan"
NAKED_DIR = os.path.join(OUTPUT_DIR, "naked")
GLYCO_DIR = os.path.join(OUTPUT_DIR, "glycosylated")


def load_trajectory(results_dir, model_name):
    """Load trajectory"""
    dcd_path = os.path.join(results_dir, f"penetration_{model_name}.dcd")
    pdb_path = os.path.join(results_dir, f"penetration_{model_name}_final.pdb")
    
    if not os.path.exists(dcd_path):
        print(f"  ❌ DCD not found: {dcd_path}")
        return None
    
    print(f"\nLoading {model_name}...")
    traj = md.load(dcd_path, top=pdb_path)
    print(f"  Frames: {traj.n_frames}, Atoms: {traj.n_atoms}")
    
    return traj


def analyze_com_distance(traj, model_name):
    """Analyze ligand COM distance to binding site"""
    print(f"\n{model_name} COM Distance Analysis:")
    
    # Select ligand and protein
    ligand_idx = traj.topology.select('not protein and not water and not (name Na or name Cl)')
    protein_ca = traj.topology.select('protein and name CA')
    
    # Find binding site center (residues near channel entrance)
    # Use residue 45 region as reference
    binding_site_residues = []
    for residue in traj.topology.residues:
        if residue.resSeq in range(40, 55):  # Near Asn45
            for atom in residue.atoms:
                if atom.name == 'CA':
                    binding_site_residues.append(atom.index)
    
    if not binding_site_residues:
        binding_site_residues = protein_ca[:20]  # Fallback
    
    print(f"  Ligand atoms: {len(ligand_idx)}")
    print(f"  Binding site reference atoms: {len(binding_site_residues)}")
    
    # Calculate COM distances
    com_distances = []
    
    for frame_idx in range(traj.n_frames):
        # Ligand COM
        lig_pos = traj.xyz[frame_idx, ligand_idx, :]
        lig_com = np.mean(lig_pos, axis=0)
        
        # Binding site COM
        bs_pos = traj.xyz[frame_idx, binding_site_residues, :]
        bs_com = np.mean(bs_pos, axis=0)
        
        dist = np.linalg.norm(lig_com - bs_com) * 10  # nm to Å
        com_distances.append(dist)
    
    com_distances = np.array(com_distances)
    
    print(f"  Initial distance: {com_distances[0]:.2f} Å")
    print(f"  Final distance: {com_distances[-1]:.2f} Å")
    print(f"  Mean distance: {np.mean(com_distances):.2f} ± {np.std(com_distances):.2f} Å")
    print(f"  Min distance: {np.min(com_distances):.2f} Å")
    
    # Did drug approach?
    approach = com_distances[0] - np.min(com_distances)
    print(f"  Closest approach: {approach:.2f} Å closer than start")
    
    return {
        'distances': com_distances,
        'initial': com_distances[0],
        'final': com_distances[-1],
        'mean': np.mean(com_distances),
        'std': np.std(com_distances),
        'min': np.min(com_distances),
        'approach': approach
    }


def analyze_glycan_contacts(traj, model_name):
    """Analyze contacts with glycan residues"""
    print(f"\n{model_name} Glycan Contact Analysis:")
    
    ligand_idx = traj.topology.select('not protein and not water and not (name Na or name Cl)')
    
    # Find glycan residues (NAG, MAN, BMA, FUC, GAL)
    glycan_idx = []
    glycan_residues = ['NAG', 'MAN', 'BMA', 'FUC', 'GAL']
    
    for atom in traj.topology.atoms():
        if atom.residue.name in glycan_residues:
            glycan_idx.append(atom.index)
    
    if not glycan_idx:
        print(f"  No glycan residues found (naked model)")
        return {'contacts': np.zeros(traj.n_frames), 'mean': 0, 'std': 0}
    
    print(f"  Glycan atoms: {len(glycan_idx)}")
    
    # Calculate contacts per frame
    contacts_per_frame = []
    cutoff = 0.5  # 5 Å
    
    for frame_idx in range(0, traj.n_frames, 5):
        frame_pos = traj.xyz[frame_idx]
        lig_pos = frame_pos[ligand_idx]
        glycan_pos = frame_pos[glycan_idx]
        
        count = 0
        for lp in lig_pos:
            dists = np.linalg.norm(glycan_pos - lp, axis=1)
            count += np.sum(dists < cutoff)
        
        contacts_per_frame.append(count)
    
    contacts = np.array(contacts_per_frame)
    
    print(f"  Mean glycan contacts: {np.mean(contacts):.1f} ± {np.std(contacts):.1f}")
    
    return {
        'contacts': contacts,
        'mean': np.mean(contacts),
        'std': np.std(contacts)
    }


def analyze_protein_contacts(traj, model_name, cutoff=0.5):
    """Analyze contacts with protein (excluding glycan)"""
    print(f"\n{model_name} Protein Contact Analysis:")
    
    ligand_idx = traj.topology.select('not protein and not water and not (name Na or name Cl)')
    protein_idx = traj.topology.select('protein and element != H')
    
    # Exclude glycan from ligand selection
    glycan_residues = ['NAG', 'MAN', 'BMA', 'FUC', 'GAL']
    ligand_idx = [i for i in ligand_idx 
                  if traj.topology.atom(i).residue.name not in glycan_residues]
    
    print(f"  Ligand atoms: {len(ligand_idx)}")
    print(f"  Protein heavy atoms: {len(protein_idx)}")
    
    contacts_per_frame = []
    
    for frame_idx in range(0, traj.n_frames, 5):
        frame_pos = traj.xyz[frame_idx]
        lig_pos = frame_pos[ligand_idx]
        prot_pos = frame_pos[protein_idx]
        
        count = 0
        for lp in lig_pos:
            dists = np.linalg.norm(prot_pos - lp, axis=1)
            count += np.sum(dists < cutoff)
        
        contacts_per_frame.append(count)
    
    contacts = np.array(contacts_per_frame)
    
    print(f"  Mean protein contacts: {np.mean(contacts):.1f} ± {np.std(contacts):.1f}")
    
    return {
        'contacts': contacts,
        'mean': np.mean(contacts),
        'std': np.std(contacts)
    }


def analyze_sasa(traj, model_name):
    """Analyze Solvent Accessible Surface Area of binding site"""
    print(f"\n{model_name} SASA Analysis:")
    
    # Find binding site residues
    binding_site_idx = []
    for residue in traj.topology.residues:
        if residue.resSeq in range(40, 60):
            for atom in residue.atoms:
                binding_site_idx.append(atom.index)
    
    if not binding_site_idx:
        print("  Could not identify binding site")
        return None
    
    # Calculate SASA for binding site
    sasa_values = md.shrake_rupley(traj, mode='residue')
    
    # Sum SASA for binding site residues
    binding_site_sasa = []
    for frame_idx in range(traj.n_frames):
        total_sasa = 0
        for residue in traj.topology.residues:
            if residue.resSeq in range(40, 60):
                total_sasa += sasa_values[frame_idx, residue.index]
        binding_site_sasa.append(total_sasa * 100)  # nm² to Å²
    
    binding_site_sasa = np.array(binding_site_sasa)
    
    print(f"  Mean binding site SASA: {np.mean(binding_site_sasa):.1f} ± {np.std(binding_site_sasa):.1f} Å²")
    
    return {
        'sasa': binding_site_sasa,
        'mean': np.mean(binding_site_sasa),
        'std': np.std(binding_site_sasa)
    }


def plot_comparison(results_naked, results_glyco):
    """Create comparison plots"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. COM Distance over time
    ax = axes[0, 0]
    time_naked = np.arange(len(results_naked['com']['distances'])) * 0.01
    time_glyco = np.arange(len(results_glyco['com']['distances'])) * 0.01
    
    ax.plot(time_naked, results_naked['com']['distances'], 'b-', alpha=0.7, 
            label=f"Naked (cancer): {results_naked['com']['mean']:.1f}Å")
    ax.plot(time_glyco, results_glyco['com']['distances'], 'r-', alpha=0.7,
            label=f"Glycosylated (normal): {results_glyco['com']['mean']:.1f}Å")
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Distance to Binding Site (Å)')
    ax.set_title('Drug Approach to GLUT1')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 2. Protein Contacts
    ax = axes[0, 1]
    time_naked = np.arange(len(results_naked['protein_contacts']['contacts'])) * 0.05
    time_glyco = np.arange(len(results_glyco['protein_contacts']['contacts'])) * 0.05
    
    ax.plot(time_naked, results_naked['protein_contacts']['contacts'], 'b-', alpha=0.7,
            label=f"Naked: {results_naked['protein_contacts']['mean']:.0f}")
    ax.plot(time_glyco, results_glyco['protein_contacts']['contacts'], 'r-', alpha=0.7,
            label=f"Glycosylated: {results_glyco['protein_contacts']['mean']:.0f}")
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Protein Contacts')
    ax.set_title('Drug-Protein Contacts')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 3. Glycan Contacts (glycosylated only)
    ax = axes[1, 0]
    if results_glyco['glycan_contacts']['mean'] > 0:
        time_glyco = np.arange(len(results_glyco['glycan_contacts']['contacts'])) * 0.05
        ax.plot(time_glyco, results_glyco['glycan_contacts']['contacts'], 'r-', alpha=0.7)
        ax.axhline(y=results_glyco['glycan_contacts']['mean'], color='r', linestyle='--', 
                   label=f"Mean: {results_glyco['glycan_contacts']['mean']:.1f}")
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Glycan Contacts')
    ax.set_title('Drug-Glycan Contacts (Normal Cell Model)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 4. Summary Bar Chart
    ax = axes[1, 1]
    metrics = ['Closest\nApproach (Å)', 'Protein\nContacts', 'Mean\nDistance (Å)']
    naked_vals = [
        results_naked['com']['approach'],
        results_naked['protein_contacts']['mean'],
        results_naked['com']['mean']
    ]
    glyco_vals = [
        results_glyco['com']['approach'],
        results_glyco['protein_contacts']['mean'],
        results_glyco['com']['mean']
    ]
    
    x = np.arange(len(metrics))
    width = 0.35
    
    ax.bar(x - width/2, naked_vals, width, label='Naked (Cancer)', color='blue', alpha=0.7)
    ax.bar(x + width/2, glyco_vals, width, label='Glycosylated (Normal)', color='red', alpha=0.7)
    
    ax.set_ylabel('Value')
    ax.set_title('Penetration Summary')
    ax.set_xticks(x)
    ax.set_xticklabels(metrics)
    ax.legend()
    
    # Add value labels
    for i, (nv, gv) in enumerate(zip(naked_vals, glyco_vals)):
        ax.text(i - width/2, nv + 1, f'{nv:.1f}', ha='center', fontsize=9)
        ax.text(i + width/2, gv + 1, f'{gv:.1f}', ha='center', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'penetration_comparison.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n📊 Plot saved: {OUTPUT_DIR}/penetration_comparison.png")


def main():
    print("=" * 70)
    print("Phase 4: Glycan Penetration Analysis")
    print("=" * 70)
    
    # Load trajectories
    traj_naked = load_trajectory(NAKED_DIR, "naked")
    traj_glyco = load_trajectory(GLYCO_DIR, "glycosylated")
    
    if traj_naked is None or traj_glyco is None:
        print("\n❌ Trajectories not found. Run simulation first.")
        return
    
    # Analyze naked model
    print("\n" + "=" * 50)
    print("ANALYZING NAKED GLUT1 (Cancer Cell)")
    print("=" * 50)
    
    results_naked = {
        'com': analyze_com_distance(traj_naked, "Naked"),
        'glycan_contacts': analyze_glycan_contacts(traj_naked, "Naked"),
        'protein_contacts': analyze_protein_contacts(traj_naked, "Naked")
    }
    
    # Analyze glycosylated model
    print("\n" + "=" * 50)
    print("ANALYZING GLYCOSYLATED GLUT1 (Normal Cell)")
    print("=" * 50)
    
    results_glyco = {
        'com': analyze_com_distance(traj_glyco, "Glycosylated"),
        'glycan_contacts': analyze_glycan_contacts(traj_glyco, "Glycosylated"),
        'protein_contacts': analyze_protein_contacts(traj_glyco, "Glycosylated")
    }
    
    # Plot comparison
    plot_comparison(results_naked, results_glyco)
    
    # Summary
    print("\n" + "=" * 70)
    print("PENETRATION TEST SUMMARY")
    print("=" * 70)
    
    print("\n┌─────────────────────────┬─────────────────┬─────────────────┐")
    print("│ Metric                  │ Naked (Cancer)  │ Glyco (Normal)  │")
    print("├─────────────────────────┼─────────────────┼─────────────────┤")
    print(f"│ Initial Distance (Å)    │ {results_naked['com']['initial']:15.1f} │ {results_glyco['com']['initial']:15.1f} │")
    print(f"│ Final Distance (Å)      │ {results_naked['com']['final']:15.1f} │ {results_glyco['com']['final']:15.1f} │")
    print(f"│ Min Distance (Å)        │ {results_naked['com']['min']:15.1f} │ {results_glyco['com']['min']:15.1f} │")
    print(f"│ Closest Approach (Å)    │ {results_naked['com']['approach']:15.1f} │ {results_glyco['com']['approach']:15.1f} │")
    print(f"│ Protein Contacts        │ {results_naked['protein_contacts']['mean']:15.1f} │ {results_glyco['protein_contacts']['mean']:15.1f} │")
    print(f"│ Glycan Contacts         │ {results_naked['glycan_contacts']['mean']:15.1f} │ {results_glyco['glycan_contacts']['mean']:15.1f} │")
    print("└─────────────────────────┴─────────────────┴─────────────────┘")
    
    # Interpretation
    print("\n📋 INTERPRETATION:")
    print("-" * 50)
    
    approach_diff = results_naked['com']['approach'] - results_glyco['com']['approach']
    if approach_diff > 5:
        print(f"✅ Drug approached {approach_diff:.1f}Å CLOSER in naked model")
        print("   → Glycan layer BLOCKS drug penetration")
    elif approach_diff > 0:
        print(f"⚠️  Drug approached {approach_diff:.1f}Å closer in naked model")
        print("   → Glycan provides some steric hindrance")
    else:
        print(f"❌ Drug approached similarly in both models")
    
    if results_glyco['glycan_contacts']['mean'] > 10:
        print(f"✅ Drug formed {results_glyco['glycan_contacts']['mean']:.0f} contacts with glycan")
        print("   → Drug gets TRAPPED in glycan layer")
    
    contact_diff = results_naked['protein_contacts']['mean'] - results_glyco['protein_contacts']['mean']
    if contact_diff > 50:
        print(f"✅ Drug made {contact_diff:.0f} MORE protein contacts in naked model")
        print("   → Better binding in cancer cells (no glycan barrier)")
    
    # Save results
    summary = {
        'naked': {
            'com_initial': float(results_naked['com']['initial']),
            'com_final': float(results_naked['com']['final']),
            'com_min': float(results_naked['com']['min']),
            'approach': float(results_naked['com']['approach']),
            'protein_contacts': float(results_naked['protein_contacts']['mean']),
            'glycan_contacts': float(results_naked['glycan_contacts']['mean'])
        },
        'glycosylated': {
            'com_initial': float(results_glyco['com']['initial']),
            'com_final': float(results_glyco['com']['final']),
            'com_min': float(results_glyco['com']['min']),
            'approach': float(results_glyco['com']['approach']),
            'protein_contacts': float(results_glyco['protein_contacts']['mean']),
            'glycan_contacts': float(results_glyco['glycan_contacts']['mean'])
        },
        'conclusion': {
            'approach_difference': float(approach_diff),
            'glycan_blocks_penetration': approach_diff > 5
        }
    }
    
    with open(os.path.join(OUTPUT_DIR, 'penetration_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n💾 Results saved: {OUTPUT_DIR}/penetration_summary.json")


if __name__ == "__main__":
    main()
