#!/usr/bin/env python
"""
MD Simulation for Model A (Basic TRIS-PEG2)
============================================
Run 1ns MD simulation using Gnina-docked pose
Compare with Model B (Arg-TRIS-PEG2)
"""

import os
import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from openmm import Platform, LangevinMiddleIntegrator, MonteCarloBarostat
from openmm.app import (
    ForceField, PDBFile, Simulation, PME, HBonds, 
    Modeller, StateDataReporter, DCDReporter, CheckpointReporter
)
from openmm.unit import kelvin, picosecond, picoseconds, femtoseconds, nanometer, atmospheres

from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator
from pdbfixer import PDBFixer

import warnings
warnings.filterwarnings('ignore')

# Paths
RECEPTOR_PDB = "/home/pjho3/projects/Drug/raw_data/4PYP.pdb"
GNINA_RESULTS = "/home/pjho3/projects/Drug/results/gnina_docking"
OUTPUT_DIR = "/home/pjho3/projects/Drug/results/phase2_gnina_md/model_a"

os.makedirs(OUTPUT_DIR, exist_ok=True)


def prepare_receptor(pdb_path):
    """Clean and prepare receptor using PDBFixer"""
    print("Preparing receptor...")
    
    # Check for cached cleaned receptor
    cached_path = "/home/pjho3/projects/Drug/results/phase2_gnina_md/cleaned_receptor.pdb"
    if os.path.exists(cached_path):
        print(f"  Using cached receptor: {cached_path}")
        return cached_path
    
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    
    output_path = cached_path
    with open(output_path, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    
    print(f"  Saved: {output_path}")
    return output_path


def get_best_pose_from_gnina(sdf_path):
    """Extract best pose from Gnina docking results"""
    print(f"Loading docked pose: {sdf_path}")
    
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
    mols = [m for m in suppl if m is not None]
    
    if not mols:
        raise ValueError(f"No valid molecules in {sdf_path}")
    
    # First pose is best
    best_mol = mols[0]
    
    # Get score if available
    if best_mol.HasProp('CNNscore'):
        score = best_mol.GetProp('CNNscore')
        print(f"  CNN score: {score}")
    
    return best_mol


def run_md_simulation(receptor_pdb, ligand_mol, output_prefix, total_ns=1.0, equil_ns=0.1):
    """Run MD simulation with OpenMM"""
    
    print("\n" + "=" * 60)
    print(f"MD Simulation: Model A (Basic TRIS-PEG2)")
    print("=" * 60)
    
    dt_fs = 2.0
    save_ps = 10.0
    
    total_steps = int(total_ns * 1e6 / dt_fs)
    equil_steps = int(equil_ns * 1e6 / dt_fs)
    save_freq = int(save_ps * 1000 / dt_fs)
    
    print(f"\nSimulation parameters:")
    print(f"  Total time: {total_ns} ns ({total_steps} steps)")
    print(f"  Equilibration: {equil_ns} ns ({equil_steps} steps)")
    print(f"  Timestep: {dt_fs} fs")
    print(f"  Save frequency: {save_ps} ps")
    
    # Save ligand as PDB
    lig_pdb_path = os.path.join(OUTPUT_DIR, "temp_lig.pdb")
    Chem.MolToPDBFile(ligand_mol, lig_pdb_path)
    print(f"\n  Ligand saved: {lig_pdb_path}")
    
    # Setup force field with GAFF for ligand
    print("\n  Setting up force field...")
    
    ligand_off = Molecule.from_rdkit(ligand_mol, allow_undefined_stereo=True)
    ligand_off.assign_partial_charges(partial_charge_method="gasteiger")
    
    gaff = GAFFTemplateGenerator(molecules=[ligand_off], forcefield='gaff-2.11')
    
    # Cache GAFF parameters
    gaff_cache = os.path.join(OUTPUT_DIR, "gaff_cache.json")
    if os.path.exists(gaff_cache):
        gaff.cache = gaff_cache
    gaff.cache = gaff_cache
    
    forcefield = ForceField(
        "amber14-all.xml",
        "amber14/tip3pfb.xml"
    )
    forcefield.registerTemplateGenerator(gaff.generator)
    
    # Pre-generate template for ligand
    print("  Generating GAFF parameters for ligand...")
    from openmm.app import PDBFile as PDBFileApp
    lig_pdb_obj = PDBFileApp(lig_pdb_path)
    try:
        _ = forcefield.createSystem(lig_pdb_obj.topology)
        print("  ✅ Ligand template generated successfully")
    except Exception as e:
        print(f"  ⚠️ Template generation warning: {e}")
    
    # Load structures
    print("  Loading receptor...")
    receptor = PDBFile(receptor_pdb)
    
    print("  Loading ligand...")
    ligand = PDBFile(lig_pdb_path)
    
    # Merge
    print("  Merging complex...")
    modeller = Modeller(receptor.topology, receptor.positions)
    modeller.add(ligand.topology, ligand.positions)
    
    # Save merged complex
    complex_pdb = os.path.join(OUTPUT_DIR, "complex_merged.pdb")
    with open(complex_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
    print(f"  Complex saved: {complex_pdb}")
    
    # Add solvent
    print("  Adding solvent (this may take a while)...")
    modeller.addSolvent(
        forcefield,
        model='tip3p',
        padding=1.0 * nanometer,
        ionicStrength=0.15,
        positiveIon='Na+',
        negativeIon='Cl-'
    )
    
    print(f"  Total atoms after solvation: {modeller.topology.getNumAtoms()}")
    
    # Create system
    print("  Creating system...")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * nanometer,
        constraints=HBonds
    )
    
    # Add barostat
    system.addForce(MonteCarloBarostat(1 * atmospheres, 310 * kelvin))
    
    # Integrator
    integrator = LangevinMiddleIntegrator(
        310 * kelvin,
        1.0 / picosecond,
        dt_fs * femtoseconds
    )
    
    # Platform
    try:
        platform = Platform.getPlatformByName("CUDA")
        properties = {'Precision': 'mixed'}
        print("  Using CUDA platform")
    except:
        platform = Platform.getPlatformByName("CPU")
        properties = {}
        print("  Using CPU platform")
    
    # Simulation
    simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)
    
    # Minimize
    print("\n  Minimizing energy...")
    simulation.minimizeEnergy(maxIterations=1000)
    
    # Equilibration
    print(f"  Equilibrating for {equil_ns} ns...")
    simulation.step(equil_steps)
    
    # Production reporters
    dcd_path = os.path.join(OUTPUT_DIR, f"{output_prefix}.dcd")
    log_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_log.csv")
    chk_path = os.path.join(OUTPUT_DIR, f"{output_prefix}.chk")
    
    simulation.reporters.append(DCDReporter(dcd_path, save_freq))
    simulation.reporters.append(StateDataReporter(
        log_path, save_freq,
        step=True, potentialEnergy=True, temperature=True,
        volume=True, speed=True
    ))
    simulation.reporters.append(CheckpointReporter(chk_path, save_freq * 10))
    
    # Production
    prod_steps = total_steps - equil_steps
    print(f"\n  Running production for {total_ns - equil_ns} ns...")
    print(f"  Progress will be saved to: {log_path}")
    
    simulation.step(prod_steps)
    
    # Save final structure
    final_pdb = os.path.join(OUTPUT_DIR, f"{output_prefix}_final.pdb")
    state = simulation.context.getState(getPositions=True)
    with open(final_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, state.getPositions(), f)
    
    print(f"\n✅ Simulation complete!")
    print(f"  Trajectory: {dcd_path}")
    print(f"  Final structure: {final_pdb}")
    
    return dcd_path, final_pdb


def main():
    print("=" * 60)
    print("Model A MD Simulation (Basic TRIS-PEG2)")
    print("=" * 60)
    
    # Prepare receptor
    receptor_pdb = prepare_receptor(RECEPTOR_PDB)
    
    # Get Model A docked pose
    model_a_sdf = os.path.join(GNINA_RESULTS, "model_a_tris_peg2_docked.sdf")
    
    if not os.path.exists(model_a_sdf):
        print(f"❌ Model A docking results not found: {model_a_sdf}")
        print("Please run Gnina docking for Model A first.")
        sys.exit(1)
    
    ligand_mol = get_best_pose_from_gnina(model_a_sdf)
    
    # Run simulation
    dcd_path, final_pdb = run_md_simulation(
        receptor_pdb,
        ligand_mol,
        output_prefix="prod_model_a",
        total_ns=1.0,
        equil_ns=0.1
    )
    
    print("\n" + "=" * 60)
    print("Model A Simulation Complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
