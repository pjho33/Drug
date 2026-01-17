# scripts/07_analyze_contacts.py
"""
Analyze ligand-pocket contacts for Phase 2 results.
Measures how often the ligand stays in contact with key binding site residues.
"""
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import os

# GLUT1 binding pocket residues (from literature and 4PYP structure)
# These are key residues that interact with glucose in the binding site
POCKET_RESIDUES = [
    "GLN161", "GLN282", "GLN283", "ASN288", "ASN317",  # Polar contacts
    "PHE26", "PHE291", "PHE379", "TRP388",  # Aromatic (cation-pi potential)
    "TRP65", "TRP412",  # Additional aromatics
    "GLU380", "GLU209",  # Charged
]

def get_pocket_indices(topology):
    """Get atom indices for pocket residues."""
    pocket_indices = []
    pocket_res_found = []
    
    for res in topology.residues:
        res_id = f"{res.name}{res.resSeq}"
        if res_id in POCKET_RESIDUES or res.name in ["PHE", "TRP", "GLN", "ASN", "GLU"]:
            # Check if this residue is near the binding site (rough filter)
            if res.resSeq in range(20, 420):  # GLUT1 TM domain range
                for atom in res.atoms:
                    pocket_indices.append(atom.index)
                if res_id not in pocket_res_found:
                    pocket_res_found.append(res_id)
    
    return pocket_indices, pocket_res_found

def analyze_contacts(output_dir="results/phase2_rep1", contact_cutoff=0.4):
    """
    Analyze ligand-pocket contacts.
    contact_cutoff: distance in nm (0.4 nm = 4 Angstrom)
    """
    results = []
    ligands = ["tripod", "glucose", "bng"]
    
    for name in ligands:
        dcd = os.path.join(output_dir, f"prod_{name}_rep1.dcd")
        pdb = os.path.join(output_dir, f"prod_{name}_rep1_final.pdb")
        
        if not os.path.exists(dcd) or not os.path.exists(pdb):
            print(f"Skipping {name}: files not found")
            continue
        
        print(f"\nAnalyzing {name}...")
        t = md.load(dcd, top=pdb)
        
        # Find ligand
        ligand_candidates = [r for r in t.topology.residues 
                            if r.name not in ("HOH", "WAT", "NA", "CL", "K", "MG")]
        ligand = ligand_candidates[-1]
        ligand_indices = [atom.index for atom in ligand.atoms]
        
        # Get pocket residues
        pocket_indices, pocket_res_found = get_pocket_indices(t.topology)
        print(f"  Found {len(pocket_res_found)} pocket residue types")
        
        # Calculate contacts per frame
        n_frames = t.n_frames
        contact_counts = []
        
        # Use mdtraj's compute_contacts for efficiency
        # Create pairs between ligand and pocket atoms
        pairs = []
        for lig_idx in ligand_indices:
            for pock_idx in pocket_indices:
                pairs.append((lig_idx, pock_idx))
        
        if len(pairs) == 0:
            print(f"  Warning: No pocket atoms found for {name}")
            continue
        
        pairs = np.array(pairs)
        
        # Compute distances in batches to save memory
        batch_size = 1000
        contacts_per_frame = np.zeros(n_frames)
        
        for i in range(0, n_frames, batch_size):
            end_idx = min(i + batch_size, n_frames)
            t_batch = t[i:end_idx]
            distances = md.compute_distances(t_batch, pairs)
            # Count contacts (distance < cutoff)
            contacts_per_frame[i:end_idx] = np.sum(distances < contact_cutoff, axis=1)
        
        avg_contacts = np.mean(contacts_per_frame)
        std_contacts = np.std(contacts_per_frame)
        
        # Calculate occupancy (% of frames with at least N contacts)
        min_contacts_threshold = 5
        occupancy = np.sum(contacts_per_frame >= min_contacts_threshold) / n_frames * 100
        
        results.append({
            "name": name,
            "avg_contacts": avg_contacts,
            "std_contacts": std_contacts,
            "occupancy": occupancy,
            "contacts_per_frame": contacts_per_frame
        })
        
        print(f"  Avg contacts: {avg_contacts:.1f} ± {std_contacts:.1f}")
        print(f"  Pocket occupancy (≥{min_contacts_threshold} contacts): {occupancy:.1f}%")
    
    # Print summary
    print("\n" + "=" * 70)
    print("Phase 2 rep1 Ligand-Pocket Contact Analysis")
    print("=" * 70)
    print(f"{'Ligand':<12} {'Avg Contacts':<15} {'Std':<10} {'Occupancy (%)':<15}")
    print("-" * 70)
    for r in results:
        print(f"{r['name']:<12} {r['avg_contacts']:<15.1f} {r['std_contacts']:<10.1f} {r['occupancy']:<15.1f}")
    print("=" * 70)
    
    # Plot contacts over time
    plt.figure(figsize=(12, 6))
    colors = {'tripod': 'royalblue', 'glucose': 'green', 'bng': 'orange'}
    
    for r in results:
        # Subsample for plotting
        contacts_sub = r['contacts_per_frame'][::10]
        time_ns = np.arange(len(contacts_sub)) * 0.1
        plt.plot(time_ns, contacts_sub, label=f"{r['name']} (avg={r['avg_contacts']:.1f})", 
                color=colors.get(r['name'], 'gray'), alpha=0.7, linewidth=1)
    
    plt.xlabel("Time (ns)")
    plt.ylabel("Number of Contacts (< 4 Å)")
    plt.title("Phase 2 rep1: Ligand-Pocket Contacts Over Time")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    plot_path = os.path.join(output_dir, "contacts_comparison.png")
    plt.savefig(plot_path, dpi=150)
    print(f"\nPlot saved to: {plot_path}")
    
    return results

if __name__ == "__main__":
    analyze_contacts()
