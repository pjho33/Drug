# scripts/08_calc_mmgbsa_phase2.py
"""
Calculate MM-GBSA binding energy for Phase 2 results.
Uses DCD trajectory + final PDB topology.
"""
import sys
import os
import numpy as np
import mdtraj as md
from openmm import *
from openmm.app import *
from openmm.unit import *
from rdkit import Chem
from openmmforcefields.generators import GAFFTemplateGenerator
from openff.toolkit.topology import Molecule

def load_ligand_for_gaff(ligand_sdf):
    """Load ligand and prepare for GAFF parameterization."""
    suppl = Chem.SDMolSupplier(ligand_sdf, removeHs=False)
    rdkit_mol = next(suppl)
    
    if rdkit_mol is None:
        raise ValueError(f"Could not load ligand from {ligand_sdf}")
    
    # Add hydrogens if needed
    if rdkit_mol.GetNumAtoms() == sum([a.GetAtomicNum() != 1 for a in rdkit_mol.GetAtoms()]):
        rdkit_mol = Chem.AddHs(rdkit_mol, addCoords=True)
    
    # Remove stereochemistry to avoid errors
    for atom in rdkit_mol.GetAtoms():
        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    for bond in rdkit_mol.GetBonds():
        if bond.GetStereo() != Chem.BondStereo.STEREONONE:
            bond.SetStereo(Chem.BondStereo.STEREONONE)
    
    off_mol = Molecule.from_rdkit(rdkit_mol, allow_undefined_stereo=True)
    return off_mol

def calculate_mmgbsa(output_dir="results/phase2_rep1", n_frames=100):
    """Calculate MM-GBSA for all ligands in Phase 2."""
    
    ligands = ["tripod", "glucose", "bng"]
    results = []
    
    for name in ligands:
        dcd = os.path.join(output_dir, f"prod_{name}_rep1.dcd")
        pdb = os.path.join(output_dir, f"prod_{name}_rep1_final.pdb")
        
        # Find ligand SDF from docking
        dock_dir = os.path.join(output_dir, f"dock_{name}", "complex_0")
        ligand_sdf = None
        for f in os.listdir(dock_dir):
            if f.startswith("rank1") and f.endswith(".sdf"):
                ligand_sdf = os.path.join(dock_dir, f)
                break
        
        if not all([os.path.exists(dcd), os.path.exists(pdb), ligand_sdf]):
            print(f"Skipping {name}: files not found")
            continue
        
        print(f"\n{'='*60}")
        print(f"Calculating MM-GBSA for {name}")
        print(f"{'='*60}")
        
        # Load trajectory with mdtraj
        print("Loading trajectory...")
        traj = md.load(dcd, top=pdb)
        total_frames = traj.n_frames
        stride = max(1, total_frames // n_frames)
        
        print(f"  Total frames: {total_frames}, using stride={stride}")
        
        # Load ligand for GAFF
        print("Setting up GAFF for ligand...")
        off_mol = load_ligand_for_gaff(ligand_sdf)
        gaff = GAFFTemplateGenerator(molecules=[off_mol])
        
        # Setup ForceField
        forcefield = ForceField('amber14-all.xml', 'implicit/obc2.xml')
        forcefield.registerTemplateGenerator(gaff.generator)
        
        # Load PDB for OpenMM
        pdb_file = PDBFile(pdb)
        
        # Create modeller and strip water/ions
        modeller = Modeller(pdb_file.topology, pdb_file.positions)
        modeller.deleteWater()
        ions = [r for r in modeller.topology.residues() if r.name in ['NA', 'CL', 'K', 'MG']]
        if ions:
            modeller.delete(ions)
        
        # Find ligand residue (last non-standard residue)
        all_residues = list(modeller.topology.residues())
        ligand_res = all_residues[-1]
        print(f"  Ligand residue: {ligand_res.name}")
        
        # Get atom indices
        ligand_atoms = [a for a in modeller.topology.atoms() if a.residue == ligand_res]
        receptor_atoms = [a for a in modeller.topology.atoms() if a.residue != ligand_res]
        
        n_lig = len(ligand_atoms)
        n_rec = len(receptor_atoms)
        print(f"  Receptor atoms: {n_rec}, Ligand atoms: {n_lig}")
        
        # Create systems
        print("Creating OpenMM systems...")
        
        # Complex system
        system_complex = forcefield.createSystem(modeller.topology, 
                                                  nonbondedMethod=NoCutoff, 
                                                  constraints=HBonds)
        
        # Receptor only
        modeller_rec = Modeller(pdb_file.topology, pdb_file.positions)
        modeller_rec.deleteWater()
        ions_rec = [r for r in modeller_rec.topology.residues() if r.name in ['NA', 'CL', 'K', 'MG']]
        if ions_rec: modeller_rec.delete(ions_rec)
        lig_chain = [list(modeller_rec.topology.chains())[-1]]
        modeller_rec.delete(lig_chain)
        system_rec = forcefield.createSystem(modeller_rec.topology, 
                                              nonbondedMethod=NoCutoff, 
                                              constraints=HBonds)
        
        # Ligand only
        modeller_lig = Modeller(pdb_file.topology, pdb_file.positions)
        modeller_lig.deleteWater()
        ions_lig = [r for r in modeller_lig.topology.residues() if r.name in ['NA', 'CL', 'K', 'MG']]
        if ions_lig: modeller_lig.delete(ions_lig)
        rec_chains = list(modeller_lig.topology.chains())[:-1]
        modeller_lig.delete(rec_chains)
        system_lig = forcefield.createSystem(modeller_lig.topology, 
                                              nonbondedMethod=NoCutoff, 
                                              constraints=HBonds)
        
        # Setup contexts
        try:
            platform = Platform.getPlatformByName('CUDA')
            prop = {'DeviceIndex': '0', 'Precision': 'mixed'}
        except:
            platform = Platform.getPlatformByName('CPU')
            prop = {}
        
        int_c = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
        int_r = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
        int_l = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
        
        ctx_complex = Context(system_complex, int_c, platform, prop)
        ctx_rec = Context(system_rec, int_r, platform, prop)
        ctx_lig = Context(system_lig, int_l, platform, prop)
        
        # Get atom indices for stripping water/ions from trajectory
        traj_atoms = list(traj.topology.atoms)
        valid_indices = [i for i, a in enumerate(traj_atoms) 
                        if a.residue.name not in ('HOH', 'WAT', 'NA', 'CL', 'K', 'MG')]
        
        # Calculate energies
        print(f"Calculating energies for {n_frames} frames...")
        energies = []
        
        for frame_idx in range(0, total_frames, stride):
            # Get positions from trajectory (nm -> OpenMM units)
            positions_nm = traj.xyz[frame_idx]
            positions = [Vec3(positions_nm[i][0], positions_nm[i][1], positions_nm[i][2]) * nanometer 
                        for i in valid_indices]
            
            if len(positions) != (n_rec + n_lig):
                continue
            
            try:
                # Complex energy
                ctx_complex.setPositions(positions)
                e_complex = ctx_complex.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalories_per_mole)
                
                # Receptor energy
                ctx_rec.setPositions(positions[:n_rec])
                e_rec = ctx_rec.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalories_per_mole)
                
                # Ligand energy
                ctx_lig.setPositions(positions[n_rec:])
                e_lig = ctx_lig.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalories_per_mole)
                
                dG = e_complex - (e_rec + e_lig)
                energies.append(dG)
                
            except Exception as e:
                continue
        
        if energies:
            avg_energy = np.mean(energies)
            std_energy = np.std(energies)
            
            results.append({
                "name": name,
                "avg": avg_energy,
                "std": std_energy,
                "n_frames": len(energies)
            })
            
            print(f"\n  MM-GBSA: {avg_energy:.2f} ± {std_energy:.2f} kcal/mol")
            
            # Save to file
            output_csv = os.path.join(output_dir, f"mmgbsa_{name}_rep1.csv")
            with open(output_csv, 'w') as f:
                f.write("Frame,BindingEnergy(kcal/mol)\n")
                for idx, e in enumerate(energies):
                    f.write(f"{idx*stride},{e:.4f}\n")
                f.write(f"AVERAGE,{avg_energy:.4f}\n")
                f.write(f"STD_DEV,{std_energy:.4f}\n")
            print(f"  Saved to: {output_csv}")
    
    # Summary
    print("\n" + "=" * 70)
    print("Phase 2 rep1 MM-GBSA Summary")
    print("=" * 70)
    print(f"{'Ligand':<12} {'MM-GBSA (kcal/mol)':<25} {'Frames':<10}")
    print("-" * 70)
    for r in results:
        print(f"{r['name']:<12} {r['avg']:.2f} ± {r['std']:.2f}{'':10} {r['n_frames']:<10}")
    print("=" * 70)
    
    # Interpretation
    print("\nInterpretation (more negative = stronger binding):")
    if results:
        sorted_results = sorted(results, key=lambda x: x['avg'])
        for i, r in enumerate(sorted_results):
            rank = i + 1
            print(f"  #{rank} {r['name']}: {r['avg']:.2f} kcal/mol")
    
    return results

if __name__ == "__main__":
    calculate_mmgbsa()
