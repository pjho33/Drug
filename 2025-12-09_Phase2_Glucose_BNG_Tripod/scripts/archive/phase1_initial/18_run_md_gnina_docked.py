#!/usr/bin/env python
"""
Phase 2: MD Simulation with Gnina-Docked Ligands
=================================================
Uses best pose from Gnina docking for MD simulation
"""

import argparse
import os
import sys
import time
import json
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


def prepare_receptor(pdb_path, output_dir):
    """Prepare GLUT1 receptor using PDBFixer"""
    print("=" * 60)
    print("Step 1: Preparing Receptor")
    print("=" * 60)
    
    cleaned_path = os.path.join(output_dir, "cleaned_receptor.pdb")
    
    # Check if already prepared
    if os.path.exists(cleaned_path):
        print(f"✅ Using existing cleaned receptor: {cleaned_path}")
        return cleaned_path
    
    fixer = PDBFixer(filename=pdb_path)
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    
    with open(cleaned_path, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    
    print(f"✅ Cleaned receptor saved: {cleaned_path}")
    return cleaned_path


def get_best_pose(docked_sdf):
    """Get best pose from Gnina docking results"""
    suppl = Chem.SDMolSupplier(docked_sdf, removeHs=False)
    mols = [m for m in suppl if m is not None]
    
    if not mols:
        raise ValueError(f"No valid molecules in {docked_sdf}")
    
    # Sort by CNN score (higher is better)
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
                      total_ns=10.0, equil_ns=1.0,
                      temperature_k=310.0, pressure_bar=1.0,
                      dt_fs=2.0, save_ps=10.0,
                      platform_name="CUDA", device_index="0"):
    """Run MD simulation with OpenMM"""
    print("\n" + "=" * 60)
    print("Step 2: Running MD Simulation")
    print("=" * 60)
    
    t0 = time.time()
    
    out_dir = os.path.dirname(output_prefix)
    os.makedirs(out_dir, exist_ok=True)
    
    traj_dcd = output_prefix + ".dcd"
    final_pdb = output_prefix + "_final.pdb"
    log_csv = output_prefix + "_log.csv"
    chk_path = output_prefix + ".chk"
    
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
    
    # Create system
    print("   Creating system...")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * nanometer,
        constraints=HBonds,
    )
    system.addForce(MonteCarloBarostat(pressure_bar * bar, temperature_k * kelvin))
    
    # Setup integrator
    integrator = LangevinMiddleIntegrator(
        temperature_k * kelvin,
        1 / picosecond,
        (dt_fs / 1000.0) * picoseconds
    )
    
    # Setup platform
    try:
        platform = Platform.getPlatformByName(platform_name)
        prop = {"DeviceIndex": device_index, "Precision": "mixed"}
        print(f"   Using platform: {platform_name}")
    except Exception:
        platform = Platform.getPlatformByName("CPU")
        prop = {}
        print("   ⚠️ CUDA not available, using CPU")
    
    # Create simulation
    simulation = Simulation(modeller.topology, system, integrator, platform, prop)
    simulation.context.setPositions(modeller.positions)
    
    # Minimize with iteration limit
    print("   Minimizing energy (max 2000 iterations)...")
    simulation.minimizeEnergy(maxIterations=2000)
    
    # Check energy
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy()
    print(f"   Energy after minimization: {energy}")
    
    # Additional minimization if needed
    if energy._value > 0:
        print("   Running additional minimization...")
        simulation.minimizeEnergy(maxIterations=2000)
        state = simulation.context.getState(getEnergy=True)
        print(f"   Energy after 2nd minimization: {state.getPotentialEnergy()}")
    
    # Set velocities
    simulation.context.setVelocitiesToTemperature(temperature_k * kelvin)
    
    # Calculate steps
    eq_steps = int((equil_ns * 1_000_000.0) / dt_fs)
    prod_steps = int((total_ns * 1_000_000.0) / dt_fs)
    report_steps = max(1, int((save_ps * 1000.0) / dt_fs))
    
    print(f"   Equilibration: {equil_ns} ns ({eq_steps} steps)")
    print(f"   Production: {total_ns} ns ({prod_steps} steps)")
    print(f"   Save interval: {save_ps} ps ({report_steps} steps)")
    
    # Setup reporters
    simulation.reporters.append(DCDReporter(traj_dcd, report_steps))
    simulation.reporters.append(StateDataReporter(
        log_csv, report_steps, step=True, potentialEnergy=True,
        temperature=True, volume=True, speed=True
    ))
    simulation.reporters.append(CheckpointReporter(chk_path, report_steps * 10))
    
    # Check for checkpoint
    if os.path.exists(chk_path):
        print(f"   📂 Loading checkpoint: {chk_path}")
        simulation.loadCheckpoint(chk_path)
    else:
        # Equilibration
        if eq_steps > 0:
            print(f"   ⏳ Running equilibration...")
            simulation.step(eq_steps)
    
    # Production
    print(f"   ⏳ Running production ({total_ns} ns)...")
    simulation.step(prod_steps)
    
    # Save final state
    state = simulation.context.getState(getPositions=True)
    with open(final_pdb, "w") as f:
        PDBFile.writeFile(simulation.topology, state.getPositions(), f)
    
    # Cleanup temp files
    for f in [temp_lig_pdb, merged_pdb]:
        if os.path.exists(f):
            try:
                os.remove(f)
            except:
                pass
    
    dt = time.time() - t0
    print(f"\n✅ Simulation complete!")
    print(f"   Time: {dt/3600:.1f} hours")
    print(f"   Trajectory: {traj_dcd}")
    print(f"   Final structure: {final_pdb}")
    print(f"   Log: {log_csv}")
    
    return traj_dcd, final_pdb, log_csv


def main():
    parser = argparse.ArgumentParser(description="MD Simulation with Gnina-Docked Ligands")
    parser.add_argument("--model", choices=["a", "b", "both"], default="both",
                        help="Model to simulate: a (basic), b (arginine), or both")
    parser.add_argument("--total_ns", type=float, default=10.0,
                        help="Total simulation time in ns")
    parser.add_argument("--equil_ns", type=float, default=1.0,
                        help="Equilibration time in ns")
    parser.add_argument("--receptor", default="/home/pjho3/projects/Drug/raw_data/4PYP.pdb",
                        help="Path to receptor PDB")
    parser.add_argument("--docking_dir", default="/home/pjho3/projects/Drug/results/gnina_docking",
                        help="Directory with Gnina docking results")
    parser.add_argument("--output_dir", default="/home/pjho3/projects/Drug/results/phase2_gnina_md",
                        help="Output directory")
    parser.add_argument("--platform", default="CUDA")
    parser.add_argument("--device", default="0")
    
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Prepare receptor once
    cleaned_receptor = prepare_receptor(args.receptor, args.output_dir)
    
    # Define models to run
    models = []
    if args.model in ["a", "both"]:
        models.append(("model_a", "model_a_tris_peg2_docked.sdf"))
    if args.model in ["b", "both"]:
        models.append(("model_b", "model_b_arg_tris_peg2_docked.sdf"))
    
    # Run simulations
    results = {}
    for model_name, docked_file in models:
        print("\n" + "=" * 60)
        print(f"Running {model_name.upper()}")
        print("=" * 60)
        
        docked_sdf = os.path.join(args.docking_dir, docked_file)
        if not os.path.exists(docked_sdf):
            print(f"❌ Docked file not found: {docked_sdf}")
            continue
        
        # Get best pose
        print(f"   Loading best pose from: {docked_sdf}")
        best_mol = get_best_pose(docked_sdf)
        
        # Create model output directory
        model_dir = os.path.join(args.output_dir, model_name)
        os.makedirs(model_dir, exist_ok=True)
        
        # Run MD
        output_prefix = os.path.join(model_dir, f"prod_{model_name}")
        try:
            traj, final, log = run_md_simulation(
                receptor_pdb=cleaned_receptor,
                ligand_mol=best_mol,
                output_prefix=output_prefix,
                total_ns=args.total_ns,
                equil_ns=args.equil_ns,
                platform_name=args.platform,
                device_index=args.device
            )
            
            results[model_name] = {
                "trajectory": traj,
                "final_pdb": final,
                "log": log
            }
        except Exception as e:
            print(f"❌ Simulation failed: {e}")
            import traceback
            traceback.print_exc()
    
    # Save results summary
    summary_path = os.path.join(args.output_dir, "simulation_summary.json")
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2)
    
    print("\n" + "=" * 60)
    print("Phase 2 Simulation Complete!")
    print("=" * 60)
    print(f"Results saved to: {args.output_dir}")
    print(f"Summary: {summary_path}")
    print("\nNext steps:")
    print("  1. Run RMSD analysis")
    print("  2. Calculate MM/PBSA binding energy")
    print("  3. Measure Cation-π distances (Model B)")


if __name__ == "__main__":
    main()
