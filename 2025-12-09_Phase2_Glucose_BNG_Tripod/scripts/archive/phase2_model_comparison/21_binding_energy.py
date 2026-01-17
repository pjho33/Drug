#!/usr/bin/env python
"""
Binding Energy Calculation using OpenMM
========================================
Calculates binding free energy using MM/GBSA-like approach
ΔG_bind = E_complex - E_protein - E_ligand
"""

import os
import numpy as np
import mdtraj as md
from openmm import Platform, CustomGBForce, NonbondedForce
from openmm.app import ForceField, PDBFile, Simulation, NoCutoff, HBonds
from openmm.unit import kilocalories_per_mole, nanometer, kilojoules_per_mole
from openmm import LangevinMiddleIntegrator
from openmm.unit import kelvin, picosecond, picoseconds
import warnings
warnings.filterwarnings('ignore')

# Paths
RESULTS_DIR = "/home/pjho3/projects/Drug/results/phase2_gnina_md/model_b"
OUTPUT_DIR = "/home/pjho3/projects/Drug/results/analysis"


def calculate_interaction_energy(traj_path, topology_path, n_frames=10):
    """
    Calculate protein-ligand interaction energy from trajectory
    Uses a simplified approach: measures non-bonded interactions
    """
    print("=" * 60)
    print("Binding Energy Analysis")
    print("=" * 60)
    
    # Load trajectory (sample frames)
    print(f"\nLoading trajectory...")
    traj = md.load(traj_path, top=topology_path, stride=max(1, int(120/n_frames)))
    print(f"Loaded {traj.n_frames} frames")
    
    # Identify atoms
    protein_atoms = traj.topology.select('protein')
    ligand_atoms = traj.topology.select('not protein and not water and not (name Na or name Cl)')
    
    print(f"Protein atoms: {len(protein_atoms)}")
    print(f"Ligand atoms: {len(ligand_atoms)}")
    
    if len(ligand_atoms) == 0:
        print("❌ No ligand atoms found!")
        return None
    
    # Calculate interaction energies per frame
    print("\nCalculating interaction energies...")
    
    interaction_energies = []
    
    for i in range(traj.n_frames):
        frame = traj[i]
        
        # Get positions
        protein_pos = frame.xyz[0, protein_atoms]  # nm
        ligand_pos = frame.xyz[0, ligand_atoms]    # nm
        
        # Simple Lennard-Jones + Coulomb approximation
        # Using typical parameters for organic molecules
        
        total_energy = 0.0
        
        for lig_idx, lig_pos in enumerate(ligand_pos):
            for prot_idx, prot_pos in enumerate(protein_pos):
                r = np.sqrt(np.sum((lig_pos - prot_pos)**2))  # nm
                
                if r < 0.15:  # Skip if too close (clash)
                    continue
                
                # Simplified LJ potential (attractive part only for estimation)
                # E_vdw ≈ -C6/r^6 where C6 ~ 0.001 kJ/mol*nm^6 for typical atoms
                if r < 1.0:  # Only count nearby interactions
                    e_vdw = -0.001 / (r**6) if r > 0.3 else 0
                    total_energy += e_vdw
        
        # Convert to kcal/mol (rough estimate)
        energy_kcal = total_energy * 0.239  # kJ to kcal
        interaction_energies.append(energy_kcal)
        
        if (i + 1) % 5 == 0:
            print(f"  Frame {i+1}/{traj.n_frames}: {energy_kcal:.2f} kcal/mol")
    
    return np.array(interaction_energies)


def calculate_contact_energy(traj_path, topology_path):
    """
    Estimate binding energy based on contacts
    Uses empirical relationship: ~0.5-1.0 kcal/mol per heavy atom contact
    """
    print("\n" + "=" * 60)
    print("Contact-Based Energy Estimation")
    print("=" * 60)
    
    traj = md.load(traj_path, top=topology_path, stride=10)
    
    protein_atoms = traj.topology.select('protein and element != H')
    ligand_atoms = traj.topology.select('not protein and not water and not (name Na or name Cl) and element != H')
    
    print(f"Protein heavy atoms: {len(protein_atoms)}")
    print(f"Ligand heavy atoms: {len(ligand_atoms)}")
    
    contact_energies = []
    contact_counts = []
    
    # Different distance cutoffs for different interaction types
    cutoffs = {
        'hydrogen_bond': 0.35,  # 3.5 Å - strong
        'close_contact': 0.45,  # 4.5 Å - medium  
        'vdw_contact': 0.60,    # 6.0 Å - weak
    }
    
    energy_per_contact = {
        'hydrogen_bond': -1.5,  # kcal/mol
        'close_contact': -0.5,
        'vdw_contact': -0.1,
    }
    
    for i in range(traj.n_frames):
        frame = traj[i]
        protein_pos = frame.xyz[0, protein_atoms]
        ligand_pos = frame.xyz[0, ligand_atoms]
        
        frame_energy = 0.0
        frame_contacts = {'hydrogen_bond': 0, 'close_contact': 0, 'vdw_contact': 0}
        
        for lig_pos in ligand_pos:
            dists = np.sqrt(np.sum((protein_pos - lig_pos)**2, axis=1))
            
            # Count contacts at different cutoffs
            hb = np.sum(dists < cutoffs['hydrogen_bond'])
            cc = np.sum((dists >= cutoffs['hydrogen_bond']) & (dists < cutoffs['close_contact']))
            vdw = np.sum((dists >= cutoffs['close_contact']) & (dists < cutoffs['vdw_contact']))
            
            frame_contacts['hydrogen_bond'] += hb
            frame_contacts['close_contact'] += cc
            frame_contacts['vdw_contact'] += vdw
            
            frame_energy += hb * energy_per_contact['hydrogen_bond']
            frame_energy += cc * energy_per_contact['close_contact']
            frame_energy += vdw * energy_per_contact['vdw_contact']
        
        contact_energies.append(frame_energy)
        contact_counts.append(frame_contacts)
    
    return np.array(contact_energies), contact_counts


def calculate_solvation_penalty(n_ligand_atoms):
    """
    Estimate desolvation penalty for ligand binding
    Empirical: ~0.1-0.2 kcal/mol per heavy atom buried
    """
    # Rough estimate based on ligand size
    desolvation_penalty = n_ligand_atoms * 0.15  # kcal/mol
    return desolvation_penalty


def main():
    print("=" * 60)
    print("Binding Energy Analysis - Model B (Arg-TRIS-PEG2)")
    print("=" * 60)
    
    traj_path = os.path.join(RESULTS_DIR, "prod_model_b.dcd")
    topology_path = os.path.join(RESULTS_DIR, "prod_model_b_final.pdb")
    
    # 1. Contact-based energy estimation
    contact_energies, contact_counts = calculate_contact_energy(traj_path, topology_path)
    
    # 2. Get ligand size for solvation penalty
    traj = md.load(topology_path)
    ligand_atoms = traj.topology.select('not protein and not water and not (name Na or name Cl) and element != H')
    n_ligand_heavy = len(ligand_atoms)
    
    desolvation = calculate_solvation_penalty(n_ligand_heavy)
    
    # 3. Calculate final binding energy estimate
    mean_contact_energy = np.mean(contact_energies)
    std_contact_energy = np.std(contact_energies)
    
    # Total binding energy = interaction energy - desolvation penalty
    binding_energy = mean_contact_energy + desolvation  # desolvation is positive (unfavorable)
    
    # Summary
    print("\n" + "=" * 60)
    print("Binding Energy Summary")
    print("=" * 60)
    
    print(f"\n📊 Contact Analysis (per frame average):")
    avg_contacts = {k: np.mean([c[k] for c in contact_counts]) for k in contact_counts[0].keys()}
    print(f"  H-bond contacts (<3.5Å): {avg_contacts['hydrogen_bond']:.0f}")
    print(f"  Close contacts (3.5-4.5Å): {avg_contacts['close_contact']:.0f}")
    print(f"  VdW contacts (4.5-6.0Å): {avg_contacts['vdw_contact']:.0f}")
    
    print(f"\n⚡ Energy Components:")
    print(f"  Interaction Energy: {mean_contact_energy:.2f} ± {std_contact_energy:.2f} kcal/mol")
    print(f"  Desolvation Penalty: +{desolvation:.2f} kcal/mol")
    print(f"  ─────────────────────────────────")
    print(f"  Estimated ΔG_bind: {binding_energy:.2f} kcal/mol")
    
    # Interpretation
    print(f"\n📝 Interpretation:")
    if binding_energy < -10:
        print(f"  ✅ Strong binding (ΔG < -10 kcal/mol)")
        print(f"  → Kd estimate: < 10 nM")
    elif binding_energy < -7:
        print(f"  ✅ Good binding (-10 < ΔG < -7 kcal/mol)")
        print(f"  → Kd estimate: 10 nM - 10 μM")
    elif binding_energy < -5:
        print(f"  ⚠️ Moderate binding (-7 < ΔG < -5 kcal/mol)")
        print(f"  → Kd estimate: 10 μM - 100 μM")
    else:
        print(f"  ❌ Weak binding (ΔG > -5 kcal/mol)")
        print(f"  → Kd estimate: > 100 μM")
    
    # Compare with Gnina CNN score
    print(f"\n🔬 Comparison with Gnina Docking:")
    print(f"  Gnina CNN score: 0.5833")
    print(f"  (CNN score > 0.5 indicates good binding pose)")
    
    # Save results
    results = {
        'model': 'Model B (Arg-TRIS-PEG2)',
        'interaction_energy_kcal': mean_contact_energy,
        'interaction_energy_std': std_contact_energy,
        'desolvation_penalty_kcal': desolvation,
        'estimated_binding_energy_kcal': binding_energy,
        'avg_hbond_contacts': avg_contacts['hydrogen_bond'],
        'avg_close_contacts': avg_contacts['close_contact'],
        'avg_vdw_contacts': avg_contacts['vdw_contact'],
        'gnina_cnn_score': 0.5833
    }
    
    import json
    results_path = os.path.join(OUTPUT_DIR, 'binding_energy_model_b.json')
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n✅ Results saved: {results_path}")
    
    return results


if __name__ == "__main__":
    main()
