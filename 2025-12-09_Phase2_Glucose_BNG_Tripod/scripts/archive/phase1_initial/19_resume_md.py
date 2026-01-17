#!/usr/bin/env python
"""
Resume MD Simulation from Checkpoint
=====================================
"""

import os
import sys
import time

from openmm import Platform, MonteCarloBarostat
from openmm.app import (
    DCDReporter, ForceField, Modeller, PDBFile, Simulation,
    StateDataReporter, CheckpointReporter, PME, HBonds
)
from openmm.unit import bar, kelvin, molar, nanometer, picosecond, picoseconds
from openmm import LangevinMiddleIntegrator

from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator
from rdkit import Chem


def resume_simulation(model_dir, remaining_steps, dt_fs=2.0, save_ps=10.0):
    """Resume simulation from checkpoint"""
    
    chk_path = os.path.join(model_dir, "prod_model_b.chk")
    traj_dcd = os.path.join(model_dir, "prod_model_b.dcd")
    log_csv = os.path.join(model_dir, "prod_model_b_log.csv")
    final_pdb = os.path.join(model_dir, "prod_model_b_final.pdb")
    
    merged_pdb = os.path.join(model_dir, "complex_merged.pdb")
    cache_file = os.path.join(model_dir, "gaff_cache.json")
    temp_lig_pdb = os.path.join(model_dir, "temp_lig.pdb")
    
    print("=" * 60)
    print("Resuming MD Simulation from Checkpoint")
    print("=" * 60)
    
    # Load ligand for GAFF
    print("   Loading ligand...")
    lig_mol = Chem.MolFromPDBFile(temp_lig_pdb, removeHs=False)
    ligand_off = Molecule.from_rdkit(lig_mol, allow_undefined_stereo=True)
    ligand_off.assign_partial_charges(partial_charge_method="gasteiger")
    gaff = GAFFTemplateGenerator(molecules=[ligand_off], cache=cache_file)
    
    # Setup force field
    print("   Setting up force field...")
    forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    forcefield.registerTemplateGenerator(gaff.generator)
    
    # Load complex
    print("   Loading complex...")
    pdb = PDBFile(merged_pdb)
    modeller = Modeller(pdb.topology, pdb.positions)
    
    # Add solvent (same as before)
    print("   Adding solvent...")
    modeller.addSolvent(
        forcefield,
        padding=1.0 * nanometer,
        ionicStrength=0.15 * molar,
    )
    
    # Create system
    print("   Creating system...")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * nanometer,
        constraints=HBonds,
    )
    system.addForce(MonteCarloBarostat(1.0 * bar, 310.0 * kelvin))
    
    # Setup integrator
    integrator = LangevinMiddleIntegrator(
        310.0 * kelvin,
        1 / picosecond,
        (dt_fs / 1000.0) * picoseconds
    )
    
    # Setup platform
    platform = Platform.getPlatformByName("CUDA")
    prop = {"DeviceIndex": "0", "Precision": "mixed"}
    
    # Create simulation
    simulation = Simulation(modeller.topology, system, integrator, platform, prop)
    simulation.context.setPositions(modeller.positions)
    
    # Load checkpoint
    print(f"   Loading checkpoint: {chk_path}")
    simulation.loadCheckpoint(chk_path)
    
    # Get current step
    state = simulation.context.getState()
    current_step = simulation.currentStep
    print(f"   Current step: {current_step}")
    
    # Setup reporters (append mode)
    report_steps = max(1, int((save_ps * 1000.0) / dt_fs))
    simulation.reporters.append(DCDReporter(traj_dcd, report_steps, append=True))
    simulation.reporters.append(StateDataReporter(
        log_csv, report_steps, step=True, potentialEnergy=True,
        temperature=True, volume=True, speed=True, append=True
    ))
    simulation.reporters.append(CheckpointReporter(chk_path, report_steps * 10))
    
    # Run remaining steps
    print(f"   Running {remaining_steps} more steps...")
    t0 = time.time()
    simulation.step(remaining_steps)
    dt = time.time() - t0
    
    # Save final state
    state = simulation.context.getState(getPositions=True)
    with open(final_pdb, "w") as f:
        PDBFile.writeFile(simulation.topology, state.getPositions(), f)
    
    print(f"\n✅ Simulation complete!")
    print(f"   Time: {dt:.1f} seconds")
    print(f"   Final structure: {final_pdb}")


if __name__ == "__main__":
    model_dir = "/home/pjho3/projects/Drug/results/phase2_gnina_md/model_b"
    
    # Calculate remaining steps
    # Total: 550,000 (50k equil + 500k prod), Current: 535,000
    # Remaining: 15,000 steps
    remaining_steps = 15000
    
    resume_simulation(model_dir, remaining_steps)
