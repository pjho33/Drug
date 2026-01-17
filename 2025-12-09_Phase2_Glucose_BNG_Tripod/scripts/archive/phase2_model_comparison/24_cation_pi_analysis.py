#!/usr/bin/env python
"""
Cation-π Interaction Analysis
==============================
Analyzes the stacking interaction between Arginine guanidinium group
and aromatic residues (TRP, PHE, TYR) in GLUT1
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


def find_arginine_atoms(traj):
    """
    Find Arginine guanidinium atoms in the ligand
    Guanidinium: NH-C(=NH2+)-NH2
    """
    ligand_idx = traj.topology.select('not protein and not water and not (name Na or name Cl)')
    
    # Find nitrogen atoms in ligand (potential guanidinium)
    ligand_nitrogens = []
    for idx in ligand_idx:
        atom = traj.topology.atom(idx)
        if atom.element.symbol == 'N':
            ligand_nitrogens.append(idx)
    
    print(f"Ligand nitrogen atoms: {len(ligand_nitrogens)}")
    for idx in ligand_nitrogens:
        atom = traj.topology.atom(idx)
        print(f"  {idx}: {atom.name}")
    
    return ligand_nitrogens


def find_aromatic_residues(traj):
    """
    Find aromatic residues (TRP, PHE, TYR) in the protein
    Returns ring center atoms for each aromatic residue
    """
    aromatic_residues = []
    
    for residue in traj.topology.residues:
        if residue.name in ['TRP', 'PHE', 'TYR']:
            # Get ring atoms
            ring_atoms = []
            for atom in residue.atoms:
                # For TRP: indole ring (CG, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3, CH2)
                # For PHE/TYR: benzene ring (CG, CD1, CD2, CE1, CE2, CZ)
                if residue.name == 'TRP':
                    if atom.name in ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']:
                        ring_atoms.append(atom.index)
                else:  # PHE, TYR
                    if atom.name in ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']:
                        ring_atoms.append(atom.index)
            
            if len(ring_atoms) > 0:
                aromatic_residues.append({
                    'residue': f"{residue.name}{residue.resSeq}",
                    'ring_atoms': ring_atoms
                })
    
    print(f"\nAromatic residues found: {len(aromatic_residues)}")
    for res in aromatic_residues[:10]:  # Show first 10
        print(f"  {res['residue']}: {len(res['ring_atoms'])} ring atoms")
    
    return aromatic_residues


def calculate_ring_center_and_normal(positions, ring_atoms):
    """
    Calculate the center and normal vector of an aromatic ring
    """
    ring_pos = positions[ring_atoms]
    center = np.mean(ring_pos, axis=0)
    
    # Calculate normal using first 3 atoms (plane fitting)
    if len(ring_atoms) >= 3:
        v1 = ring_pos[1] - ring_pos[0]
        v2 = ring_pos[2] - ring_pos[0]
        normal = np.cross(v1, v2)
        normal = normal / np.linalg.norm(normal)
    else:
        normal = np.array([0, 0, 1])
    
    return center, normal


def analyze_cation_pi(traj_path, topology_path):
    """
    Analyze cation-π interactions between Arginine and aromatic residues
    """
    print("=" * 60)
    print("Cation-π Interaction Analysis")
    print("=" * 60)
    
    # Load trajectory
    print("\nLoading trajectory...")
    traj = md.load(traj_path, top=topology_path, stride=1)
    print(f"Loaded {traj.n_frames} frames")
    
    # Find Arginine nitrogens in ligand
    ligand_nitrogens = find_arginine_atoms(traj)
    
    # Find aromatic residues
    aromatic_residues = find_aromatic_residues(traj)
    
    # Focus on TRP412 specifically (mentioned in H-bond analysis)
    trp412 = None
    for res in aromatic_residues:
        if 'TRP412' in res['residue']:
            trp412 = res
            break
    
    if trp412:
        print(f"\n🎯 Found TRP412 with {len(trp412['ring_atoms'])} ring atoms")
    
    # Calculate distances and angles for each frame
    results = {
        'distances': [],
        'angles': [],
        'stacking_quality': []
    }
    
    print("\nAnalyzing cation-π geometry...")
    
    for frame_idx in range(traj.n_frames):
        positions = traj.xyz[frame_idx]  # nm
        
        frame_results = []
        
        for aro_res in aromatic_residues:
            ring_center, ring_normal = calculate_ring_center_and_normal(
                positions, aro_res['ring_atoms']
            )
            
            for n_idx in ligand_nitrogens:
                n_pos = positions[n_idx]
                
                # Distance from N to ring center
                dist_vec = n_pos - ring_center
                distance = np.linalg.norm(dist_vec) * 10  # nm to Å
                
                # Angle between N-center vector and ring normal
                # For ideal stacking, this should be close to 0° or 180°
                cos_angle = np.abs(np.dot(dist_vec / np.linalg.norm(dist_vec), ring_normal))
                angle = np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))
                
                # Stacking quality: good if distance 3-5 Å and angle < 30°
                is_stacking = (3.0 < distance < 6.0) and (angle < 45)
                
                frame_results.append({
                    'residue': aro_res['residue'],
                    'n_atom': n_idx,
                    'distance': distance,
                    'angle': angle,
                    'is_stacking': is_stacking
                })
        
        results['distances'].append(frame_results)
    
    # Analyze TRP412 specifically
    print("\n" + "=" * 60)
    print("TRP412 - Ligand Cation-π Analysis")
    print("=" * 60)
    
    trp412_distances = []
    trp412_angles = []
    
    for frame_idx, frame_results in enumerate(results['distances']):
        for r in frame_results:
            if 'TRP412' in r['residue']:
                trp412_distances.append(r['distance'])
                trp412_angles.append(r['angle'])
    
    if trp412_distances:
        trp412_distances = np.array(trp412_distances)
        trp412_angles = np.array(trp412_angles)
        
        print(f"\nDistance to TRP412 ring center:")
        print(f"  Mean: {trp412_distances.mean():.2f} Å")
        print(f"  Min: {trp412_distances.min():.2f} Å")
        print(f"  Max: {trp412_distances.max():.2f} Å")
        
        print(f"\nAngle with TRP412 ring plane:")
        print(f"  Mean: {trp412_angles.mean():.1f}°")
        print(f"  Min: {trp412_angles.min():.1f}°")
        print(f"  Max: {trp412_angles.max():.1f}°")
        
        # Check stacking criteria
        stacking_frames = np.sum((trp412_distances < 6.0) & (trp412_angles < 45))
        total_measurements = len(trp412_distances)
        stacking_percentage = stacking_frames / total_measurements * 100
        
        print(f"\n🔬 Cation-π Stacking Analysis:")
        print(f"  Stacking criteria: distance < 6 Å, angle < 45°")
        print(f"  Stacking observed: {stacking_percentage:.1f}% of measurements")
        
        if stacking_percentage > 50:
            print(f"  ✅ Strong cation-π interaction with TRP412!")
        elif stacking_percentage > 20:
            print(f"  ⚠️ Moderate cation-π interaction with TRP412")
        else:
            print(f"  ❌ Weak or no cation-π stacking with TRP412")
    
    # Find all aromatic residues with potential cation-π
    print("\n" + "=" * 60)
    print("All Aromatic Residues - Cation-π Summary")
    print("=" * 60)
    
    residue_summary = {}
    for frame_results in results['distances']:
        for r in frame_results:
            res_name = r['residue']
            if res_name not in residue_summary:
                residue_summary[res_name] = {'distances': [], 'angles': [], 'stacking': 0}
            residue_summary[res_name]['distances'].append(r['distance'])
            residue_summary[res_name]['angles'].append(r['angle'])
            if r['is_stacking']:
                residue_summary[res_name]['stacking'] += 1
    
    # Sort by stacking count
    sorted_residues = sorted(residue_summary.items(), 
                            key=lambda x: x[1]['stacking'], reverse=True)
    
    print(f"\nTop aromatic residues by cation-π stacking:")
    print(f"{'Residue':<12} {'Mean Dist (Å)':<15} {'Mean Angle (°)':<15} {'Stacking %':<12}")
    print("-" * 55)
    
    for res_name, data in sorted_residues[:10]:
        mean_dist = np.mean(data['distances'])
        mean_angle = np.mean(data['angles'])
        stacking_pct = data['stacking'] / len(data['distances']) * 100
        print(f"{res_name:<12} {mean_dist:<15.2f} {mean_angle:<15.1f} {stacking_pct:<12.1f}")
    
    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Distance distribution for top residues
    top_residues = [r[0] for r in sorted_residues[:5]]
    for res_name in top_residues:
        data = residue_summary[res_name]
        axes[0].hist(data['distances'], bins=20, alpha=0.5, label=res_name)
    
    axes[0].axvline(x=4.0, color='r', linestyle='--', label='Ideal distance (4 Å)')
    axes[0].set_xlabel('Distance to Ring Center (Å)')
    axes[0].set_ylabel('Count')
    axes[0].set_title('Cation-π Distance Distribution')
    axes[0].legend()
    
    # Angle distribution
    for res_name in top_residues:
        data = residue_summary[res_name]
        axes[1].hist(data['angles'], bins=20, alpha=0.5, label=res_name)
    
    axes[1].axvline(x=30, color='r', linestyle='--', label='Stacking threshold (30°)')
    axes[1].set_xlabel('Angle with Ring Plane (°)')
    axes[1].set_ylabel('Count')
    axes[1].set_title('Cation-π Angle Distribution')
    axes[1].legend()
    
    plt.tight_layout()
    plot_path = os.path.join(OUTPUT_DIR, 'cation_pi_analysis.png')
    plt.savefig(plot_path, dpi=150)
    print(f"\n✅ Cation-π plot saved: {plot_path}")
    plt.close()
    
    return results, residue_summary


def extract_binding_site_pdb(topology_path, output_path):
    """
    Extract binding site region for visualization
    """
    print("\n" + "=" * 60)
    print("Extracting Binding Site for Visualization")
    print("=" * 60)
    
    traj = md.load(topology_path)
    
    # Get ligand atoms
    ligand_idx = traj.topology.select('not protein and not water and not (name Na or name Cl)')
    
    if len(ligand_idx) == 0:
        print("⚠️ No ligand found!")
        return None
    
    # Get ligand center
    ligand_center = np.mean(traj.xyz[0, ligand_idx], axis=0)
    
    # Select protein atoms within 10 Å of ligand center
    protein_idx = traj.topology.select('protein')
    protein_pos = traj.xyz[0, protein_idx]
    
    distances = np.sqrt(np.sum((protein_pos - ligand_center)**2, axis=1))
    nearby_protein = protein_idx[distances < 1.0]  # 10 Å = 1.0 nm
    
    # Combine ligand and nearby protein
    binding_site_idx = np.concatenate([nearby_protein, ligand_idx])
    
    # Extract and save
    binding_site = traj.atom_slice(binding_site_idx)
    binding_site.save_pdb(output_path)
    
    print(f"Binding site atoms: {len(binding_site_idx)}")
    print(f"  Protein atoms: {len(nearby_protein)}")
    print(f"  Ligand atoms: {len(ligand_idx)}")
    print(f"✅ Saved: {output_path}")
    
    # Also identify key residues
    key_residues = set()
    for idx in nearby_protein:
        atom = traj.topology.atom(idx)
        key_residues.add(f"{atom.residue.name}{atom.residue.resSeq}")
    
    print(f"\nKey residues in binding site ({len(key_residues)}):")
    aromatic = [r for r in key_residues if any(x in r for x in ['TRP', 'PHE', 'TYR'])]
    charged = [r for r in key_residues if any(x in r for x in ['ARG', 'LYS', 'ASP', 'GLU'])]
    
    print(f"  Aromatic: {', '.join(sorted(aromatic))}")
    print(f"  Charged: {', '.join(sorted(charged))}")
    
    return output_path


def main():
    print("=" * 60)
    print("Cation-π Analysis - Model B (Arg-TRIS-PEG2)")
    print("=" * 60)
    
    traj_path = os.path.join(RESULTS_DIR, "prod_model_b.dcd")
    topology_path = os.path.join(RESULTS_DIR, "prod_model_b_final.pdb")
    
    # 1. Analyze cation-π interactions
    results, residue_summary = analyze_cation_pi(traj_path, topology_path)
    
    # 2. Extract binding site for visualization
    binding_site_pdb = os.path.join(OUTPUT_DIR, "binding_site_model_b.pdb")
    extract_binding_site_pdb(topology_path, binding_site_pdb)
    
    # Save results
    summary = {
        'model': 'Model B (Arg-TRIS-PEG2)',
        'top_cation_pi_residues': []
    }
    
    sorted_residues = sorted(residue_summary.items(), 
                            key=lambda x: x[1]['stacking'], reverse=True)
    
    for res_name, data in sorted_residues[:5]:
        summary['top_cation_pi_residues'].append({
            'residue': res_name,
            'mean_distance': float(np.mean(data['distances'])),
            'mean_angle': float(np.mean(data['angles'])),
            'stacking_percentage': float(data['stacking'] / len(data['distances']) * 100)
        })
    
    results_path = os.path.join(OUTPUT_DIR, 'cation_pi_model_b.json')
    with open(results_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n✅ Results saved: {results_path}")
    print(f"\n📁 Files for visualization:")
    print(f"  - {binding_site_pdb}")
    print(f"  - Open in PyMOL or VMD to visualize Arg-TRP412 stacking")


if __name__ == "__main__":
    main()
