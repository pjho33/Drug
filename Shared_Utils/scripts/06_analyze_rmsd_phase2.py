# scripts/06_analyze_rmsd_phase2.py
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import os

def analyze_phase2_rmsd(output_dir="results/phase2_rep1"):
    """Analyze ligand RMSD for Phase 2 rep1 results."""
    
    results = []
    ligands = ["tripod", "glucose", "bng"]
    
    for name in ligands:
        dcd = os.path.join(output_dir, f"prod_{name}_rep1.dcd")
        pdb = os.path.join(output_dir, f"prod_{name}_rep1_final.pdb")
        
        if not os.path.exists(dcd) or not os.path.exists(pdb):
            print(f"Skipping {name}: files not found")
            continue
        
        print(f"Loading {name}...")
        t = md.load(dcd, top=pdb)
        
        # Find ligand (last non-water/ion residue)
        ligand_candidates = [r for r in t.topology.residues 
                            if r.name not in ("HOH", "WAT", "NA", "CL", "K", "MG")]
        ligand = ligand_candidates[-1]
        ligand_indices = [atom.index for atom in ligand.atoms]
        
        # Superpose on protein backbone
        prot_indices = t.topology.select("protein and backbone")
        t.superpose(t, frame=0, atom_indices=prot_indices)
        
        # Calculate ligand RMSD (nm -> Angstrom)
        rmsd = md.rmsd(t, t, frame=0, atom_indices=ligand_indices) * 10.0
        
        avg = np.mean(rmsd)
        std = np.std(rmsd)
        max_rmsd = np.max(rmsd)
        
        results.append({
            "name": name,
            "avg": avg,
            "std": std,
            "max": max_rmsd,
            "frames": len(rmsd),
            "rmsd": rmsd
        })
        print(f"  {name}: avg={avg:.2f} A, std={std:.2f} A, max={max_rmsd:.2f} A")
    
    # Print summary table
    print()
    print("=" * 60)
    print("Phase 2 rep1 Ligand RMSD Summary")
    print("=" * 60)
    print(f"{'Ligand':<12} {'Avg RMSD (A)':<15} {'Std (A)':<12} {'Max (A)':<12} {'Frames':<10}")
    print("-" * 60)
    for r in results:
        print(f"{r['name']:<12} {r['avg']:<15.2f} {r['std']:<12.2f} {r['max']:<12.2f} {r['frames']:<10}")
    print("=" * 60)
    
    # Interpretation
    print("\nInterpretation:")
    for r in results:
        if r['avg'] < 2.0:
            status = "🔒 Stable - tightly bound"
        elif r['avg'] > 5.0:
            status = "💃 Unstable - may have left pocket"
        else:
            status = "⚖️ Moderate - adjusting position"
        print(f"  {r['name']}: {status}")
    
    # Plot RMSD comparison
    plt.figure(figsize=(12, 6))
    colors = {'tripod': 'royalblue', 'glucose': 'green', 'bng': 'orange'}
    
    for r in results:
        # Subsample for plotting (every 10th frame)
        rmsd_sub = r['rmsd'][::10]
        time_ns = np.arange(len(rmsd_sub)) * 0.1  # 10ps * 10 = 100ps = 0.1ns
        plt.plot(time_ns, rmsd_sub, label=f"{r['name']} (avg={r['avg']:.2f} A)", 
                color=colors.get(r['name'], 'gray'), alpha=0.8, linewidth=1)
    
    plt.xlabel("Time (ns)")
    plt.ylabel("Ligand RMSD (Å)")
    plt.title("Phase 2 rep1: Ligand Stability Comparison (100ns)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    plot_path = os.path.join(output_dir, "rmsd_comparison.png")
    plt.savefig(plot_path, dpi=150)
    print(f"\nPlot saved to: {plot_path}")
    
    return results

if __name__ == "__main__":
    analyze_phase2_rmsd()
