#!/usr/bin/env python
"""
Resume Model A MD Simulation from checkpoint
=============================================
Resumes simulation from last checkpoint if available,
otherwise starts fresh with improved stability settings.
"""

import os
import sys
import time
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem

from openmm import Platform, MonteCarloBarostat
from openmm.app import (
    DCDReporter, ForceField, Modeller, PDBFile, Simulation,
    StateDataReporter, CheckpointReporter, PME, HBonds
)
from openmm.unit import bar, kelvin, molar, nanometer, picosecond, picoseconds
from openmm import LangevinMiddleIntegrator

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


def get_best_pose(docked_sdf):
    """Get best pose from Gnina docking results"""
    suppl = Chem.SDMolSupplier(docked_sdf, removeHs=False)
    mols = [m for m in suppl if m is not None]
    
    if not mols:
        raise ValueError(f"No valid molecules in {docked_sdf}")
    
    scored = []
    for mol in mols:
        if mol.HasProp("CNNscore"):
            score = float(mol.GetProp("CNNscore"))
            scored.append((score, mol))
    
    if scored:
        scored.sort(reverse=True)
        best_mol = scored[0][1]
        print(f"   Best pose CNN score: {scored[0][0]:.4f}")
    else:
        best_mol = mols[0]
        print("   Using first pose (no CNN scores)")
    
    return best_mol


def run_md_simulation(receptor_pdb, ligand_mol, output_prefix,
                      total_ns=1.0, equil_ns=0.1,
                      temperature_k=310.0, dt_fs=2.0, save_ps=10.0):
    """Run MD simulation with OpenMM - with checkpoint support"""
    print("\n" + "=" * 60)
    print("Running MD Simulation: Model A (Basic TRIS-PEG2)")
    print("=" * 60)
    
    t0 = time.time()
    
    out_dir = os.path.dirname(output_prefix)
    os.makedirs(out_dir, exist_ok=True)
    
    traj_dcd = output_prefix + ".dcd"
    final_pdb = output_prefix + "_final.pdb"
    log_csv = output_prefix + "_log.csv"
    chk_path = output_prefix + ".chk"
    
    # Check for existing checkpoint
    resume_from_checkpoint = os.path.exists(chk_path)
    
    # Prepare ligand
    print("   Preparing ligand...")
    ligand_rdkit = Chem.AddHs(ligand_mol, addCoords=True)
    
    # Create OpenFF molecule for GAFF
    ligand_off = Molecule.from_rdkit(ligand_rdkit, allow_undefined_stereo=True)
    ligand_off.assign_partial_charges(partial_charge_method="gasteiger")
    
    cache_file = os.path.join(out_dir, "gaff_cache.json")
    gaff = GAFFTemplateGenerator(molecules=[ligand_off], cache=cache_file)
    
    # Save ligand as PDB
    temp_lig_pdb = os.path.join(out_dir, "temp_lig.pdb")
    Chem.MolToPDBFile(ligand_rdkit, temp_lig_pdb)
    
    # Merge protein and ligand
    print("   Merging complex...")
    lines_prot = [l for l in open(receptor_pdb).readlines() 
                  if not l.startswith("END") and not l.startswith("CONECT")]
    lines_lig = [l for l in open(temp_lig_pdb).readlines() 
                 if l.startswith("HETATM") or l.startswith("ATOM") or l.startswith("CONECT")]
    
    merged_pdb = os.path.join(out_dir, "complex_merged.pdb")
    with open(merged_pdb, "w") as f:
        f.writelines(lines_prot)
        f.writelines(lines_lig)
        f.write("END\n")
    
    # Setup force field
    print("   Setting up force field...")
    forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    forcefield.registerTemplateGenerator(gaff.generator)
    
    # Load complex
    pdb = PDBFile(merged_pdb)
    modeller = Modeller(pdb.topology, pdb.positions)
    
    # Add solvent
    print("   Adding solvent (this may take a few minutes)...")
    modeller.addSolvent(
        forcefield,
        padding=1.0 * nanometer,
        ionicStrength=0.15 * molar,
    )
    
    print(f"   Total atoms: {modeller.topology.getNumAtoms()}")
    
    # Create system
    print("   Creating system...")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * nanometer,
        constraints=HBonds
    )
    
    # Add barostat
    system.addForce(MonteCarloBarostat(1 * bar, temperature_k * kelvin))
    
    # Integrator
    integrator = LangevinMiddleIntegrator(
        temperature_k * kelvin,
        1.0 / picosecond,
        dt_fs * 0.001 * picoseconds
    )
    
    # Platform
    try:
        platform = Platform.getPlatformByName("CUDA")
        properties = {'DeviceIndex': '0', 'Precision': 'mixed'}
        print("   Using CUDA platform")
    except:
        platform = Platform.getPlatformByName("CPU")
        properties = {}
        print("   Using CPU platform")
    
    # Simulation
    simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)
    
    # Calculate steps
    total_steps = int(total_ns * 1e6 / dt_fs)
    equil_steps = int(equil_ns * 1e6 / dt_fs)
    save_freq = int(save_ps * 1000 / dt_fs)
    prod_steps = total_steps - equil_steps
    
    # Resume or start fresh
    if resume_from_checkpoint:
        print(f"\n   📂 Resuming from checkpoint: {chk_path}")
        simulation.loadCheckpoint(chk_path)
        current_step = simulation.currentStep
        remaining_steps = total_steps - current_step
        print(f"   Current step: {current_step}")
        print(f"   Remaining steps: {remaining_steps}")
        
        if remaining_steps <= 0:
            print("   ✅ Simulation already complete!")
            return traj_dcd, final_pdb
    else:
        print("\n   🆕 Starting fresh simulation...")
        
        # Minimize - more aggressive to avoid NaN
        print("   Minimizing energy (phase 1: rough)...")
        simulation.minimizeEnergy(maxIterations=5000, tolerance=100)
        
        print("   Minimizing energy (phase 2: fine)...")
        simulation.minimizeEnergy(maxIterations=10000, tolerance=10)
        
        # Gradual heating to avoid NaN
        print("   Gradual heating...")
        simulation.context.setVelocitiesToTemperature(50 * kelvin)
        for temp in [50, 100, 150, 200, 250, 300, temperature_k]:
            integrator.setTemperature(temp * kelvin)
            simulation.step(1000)
            print(f"     Heated to {temp} K")
        
        print(f"\n   Simulation parameters:")
        print(f"     Total: {total_ns} ns ({total_steps} steps)")
        print(f"     Equilibration: {equil_ns} ns ({equil_steps} steps)")
        print(f"     Production: {total_ns - equil_ns} ns ({prod_steps} steps)")
        
        # Equilibration
        print(f"\n   Equilibrating for {equil_ns} ns...")
        chunk_size = 10000
        for i in range(0, equil_steps, chunk_size):
            steps_to_run = min(chunk_size, equil_steps - i)
            simulation.step(steps_to_run)
            if (i + steps_to_run) % 50000 == 0:
                print(f"     Equilibration progress: {(i + steps_to_run) / equil_steps * 100:.0f}%")
        
        remaining_steps = prod_steps
    
    # Production reporters - append mode for resume
    if resume_from_checkpoint:
        simulation.reporters.append(DCDReporter(traj_dcd, save_freq, append=True))
        simulation.reporters.append(StateDataReporter(
            log_csv, save_freq,
            step=True, potentialEnergy=True, temperature=True,
            volume=True, speed=True, append=True
        ))
    else:
        simulation.reporters.append(DCDReporter(traj_dcd, save_freq))
        simulation.reporters.append(StateDataReporter(
            log_csv, save_freq,
            step=True, potentialEnergy=True, temperature=True,
            volume=True, speed=True
        ))
    
    # Checkpoint every 5000 steps (~10 ps) for more frequent saves
    simulation.reporters.append(CheckpointReporter(chk_path, 5000))
    
    # Production - run in chunks with progress
    print(f"   Running production...")
    chunk_size = 50000  # 100 ps chunks
    steps_done = 0
    
    while steps_done < remaining_steps:
        steps_to_run = min(chunk_size, remaining_steps - steps_done)
        simulation.step(steps_to_run)
        steps_done += steps_to_run
        progress = steps_done / remaining_steps * 100
        print(f"     Production progress: {progress:.1f}%")
    
    # Save final structure
    state = simulation.context.getState(getPositions=True)
    with open(final_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, state.getPositions(), f)
    
    elapsed = time.time() - t0
    print(f"\n✅ Simulation complete!")
    print(f"   Time: {elapsed:.1f} seconds")
    print(f"   Trajectory: {traj_dcd}")
    print(f"   Final structure: {final_pdb}")
    
    return traj_dcd, final_pdb


def main():
    print("=" * 60)
    print("Model A MD Simulation (Resume/Start)")
    print("=" * 60)
    
    # Use cleaned receptor from Model B
    receptor_pdb = "/home/pjho3/projects/Drug/results/phase2_gnina_md/cleaned_receptor.pdb"
    
    if not os.path.exists(receptor_pdb):
        print("❌ Cleaned receptor not found. Run Model B simulation first.")
        sys.exit(1)
    
    # Get Model A docked pose
    model_a_sdf = os.path.join(GNINA_RESULTS, "model_a_tris_peg2_docked.sdf")
    
    if not os.path.exists(model_a_sdf):
        print(f"❌ Model A docking results not found: {model_a_sdf}")
        sys.exit(1)
    
    print(f"\nLoading docked pose: {model_a_sdf}")
    ligand_mol = get_best_pose(model_a_sdf)
    
    # Run simulation
    output_prefix = os.path.join(OUTPUT_DIR, "prod_model_a")
    
    dcd_path, final_pdb = run_md_simulation(
        receptor_pdb,
        ligand_mol,
        output_prefix,
        total_ns=1.0,
        equil_ns=0.1
    )
    
    print("\n" + "=" * 60)
    print("Model A Simulation Complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
