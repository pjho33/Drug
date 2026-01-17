#!/usr/bin/env python
"""
Phase 2: Short-Tripod MD Simulation
====================================
Compare Model A (Basic TRIS) vs Model B (Arginine-TRIS) binding to GLUT1

Steps:
1. Prepare receptor (GLUT1 4PYP)
2. Manual docking of ligands at channel entrance
3. Run 100ns MD simulation
4. Analyze RMSD and binding energy
"""

import argparse
import os
import sys
import time
import json
import numpy as np

from rdkit import Chem, Geometry
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
    
    fixer = PDBFixer(filename=pdb_path)
    
    # Remove heterogens (ligands, water)
    fixer.removeHeterogens(keepWater=False)
    
    # Find and add missing residues
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)  # pH 7.4
    
    # Save cleaned receptor
    cleaned_path = os.path.join(output_dir, "cleaned_receptor.pdb")
    with open(cleaned_path, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    
    print(f"✅ Cleaned receptor saved: {cleaned_path}")
    return cleaned_path


def get_binding_site_center(receptor_pdb):
    """Get GLUT1 binding site center from known residues"""
    from Bio.PDB import PDBParser
    
    # GLUT1 glucose binding site residues (4PYP)
    BINDING_RESIDUES = [34, 161, 282, 283, 288, 317, 388, 412]
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('receptor', receptor_pdb)
    
    binding_coords = []
    all_coords = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    coord = residue['CA'].get_coord()
                    all_coords.append(coord)
                    if residue.get_id()[1] in BINDING_RESIDUES:
                        binding_coords.append(coord)
    
    all_coords = np.array(all_coords)
    center = all_coords.mean(axis=0)
    z_max = all_coords[:, 2].max()
    
    if binding_coords:
        binding_center = np.array(binding_coords).mean(axis=0)
    else:
        binding_center = center
    
    return center, binding_center, z_max


def place_ligand_at_entrance(ligand_sdf, receptor_pdb, output_sdf):
    """Place ligand at GLUT1 channel entrance"""
    print("\n" + "=" * 60)
    print("Step 2: Placing Ligand at Channel Entrance")
    print("=" * 60)
    
    # Load ligand
    suppl = Chem.SDMolSupplier(ligand_sdf, removeHs=False)
    mol = next(suppl)
    if mol is None:
        raise ValueError(f"Failed to load ligand from {ligand_sdf}")
    
    # Add hydrogens if needed
    if mol.GetNumAtoms() == sum([1 for a in mol.GetAtoms() if a.GetAtomicNum() != 1]):
        mol = Chem.AddHs(mol, addCoords=True)
    
    conf = mol.GetConformer()
    
    # Get receptor info
    center, binding_center, z_max = get_binding_site_center(receptor_pdb)
    print(f"   Receptor center: ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})")
    print(f"   Binding site: ({binding_center[0]:.1f}, {binding_center[1]:.1f}, {binding_center[2]:.1f})")
    print(f"   Z_max: {z_max:.1f} Å")
    
    # Center ligand at origin
    lig_coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    lig_center = lig_coords.mean(axis=0)
    
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Geometry.Point3D(
            float(pos.x - lig_center[0]),
            float(pos.y - lig_center[1]),
            float(pos.z - lig_center[2])))
    
    # Find arm tips (furthest atoms from center - likely glucose oxygens)
    lig_coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    distances = np.linalg.norm(lig_coords, axis=1)
    arm_z_avg = lig_coords[distances > np.percentile(distances, 90), 2].mean()
    
    # Flip if arms point up (we want arms pointing down into channel)
    if arm_z_avg > 0:
        print("   🔄 Flipping ligand (arms down)")
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, Geometry.Point3D(float(pos.x), float(-pos.y), float(-pos.z)))
    
    # Position: closer to channel entrance for stability
    # Short-tripod has ~10Å radius, place core 5Å above z_max
    target_z = z_max + 5
    target_pos = np.array([center[0], center[1], target_z])
    
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Geometry.Point3D(
            float(pos.x + target_pos[0]),
            float(pos.y + target_pos[1]),
            float(pos.z + target_pos[2])))
    
    # Verify placement
    final_coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    print(f"   Ligand placed at Z={final_coords[:, 2].mean():.1f} Å")
    print(f"   Ligand Z range: {final_coords[:, 2].min():.1f} ~ {final_coords[:, 2].max():.1f} Å")
    
    # Save
    writer = Chem.SDWriter(output_sdf)
    writer.write(mol)
    writer.close()
    print(f"   ✅ Saved: {output_sdf}")
    
    return output_sdf


def run_md_simulation(receptor_pdb, ligand_sdf, output_prefix,
                      total_ns=100.0, equil_ns=1.0,
                      temperature_k=310.0, pressure_bar=1.0,
                      dt_fs=2.0, save_ps=10.0,
                      platform_name="CUDA", device_index="0"):
    """Run MD simulation with OpenMM"""
    print("\n" + "=" * 60)
    print("Step 3: Running MD Simulation")
    print("=" * 60)
    
    t0 = time.time()
    
    out_dir = os.path.dirname(output_prefix)
    os.makedirs(out_dir, exist_ok=True)
    
    traj_dcd = output_prefix + ".dcd"
    final_pdb = output_prefix + "_final.pdb"
    log_csv = output_prefix + "_log.csv"
    chk_path = output_prefix + ".chk"
    
    # Load and prepare ligand
    print("   Loading ligand...")
    suppl = Chem.SDMolSupplier(ligand_sdf, removeHs=False)
    ligand_rdkit = next(suppl)
    if ligand_rdkit is None:
        raise ValueError("Failed to load ligand")
    
    # Ensure 3D coords and hydrogens
    if ligand_rdkit.GetNumAtoms() == sum([1 for a in ligand_rdkit.GetAtoms() if a.GetAtomicNum() != 1]):
        ligand_rdkit = Chem.AddHs(ligand_rdkit, addCoords=True)
    
    # Create OpenFF molecule for GAFF
    ligand_off = Molecule.from_rdkit(ligand_rdkit, allow_undefined_stereo=True)
    
    # Setup GAFF template generator with gasteiger charges (faster for large molecules)
    cache_file = os.path.join(out_dir, "gaff_cache.json")
    
    # Assign gasteiger charges before GAFF (am1bcc fails for large molecules)
    ligand_off.assign_partial_charges(partial_charge_method="gasteiger")
    
    gaff = GAFFTemplateGenerator(molecules=[ligand_off], cache=cache_file)
    
    # Merge protein and ligand
    print("   Merging complex...")
    temp_lig_pdb = os.path.join(out_dir, "temp_lig.pdb")
    Chem.MolToPDBFile(ligand_rdkit, temp_lig_pdb)
    
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
    
    # Add solvent (reduced padding for faster setup)
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
    
    # Check for NaN after minimization
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    energy = state.getPotentialEnergy()
    print(f"   Energy after minimization: {energy}")
    
    # Additional minimization if energy is still high
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
    parser = argparse.ArgumentParser(description="Phase 2: Short-Tripod MD Simulation")
    parser.add_argument("--model", choices=["a", "b", "both"], default="both",
                        help="Model to simulate: a (basic), b (arginine), or both")
    parser.add_argument("--peg", type=int, choices=[1, 2], default=2,
                        help="PEG units (1 or 2)")
    parser.add_argument("--total_ns", type=float, default=100.0,
                        help="Total simulation time in ns")
    parser.add_argument("--equil_ns", type=float, default=1.0,
                        help="Equilibration time in ns")
    parser.add_argument("--receptor", default="/home/pjho3/projects/Drug/raw_data/4PYP.pdb",
                        help="Path to receptor PDB")
    parser.add_argument("--output_dir", default="/home/pjho3/projects/Drug/results/phase2_short_tripod",
                        help="Output directory")
    parser.add_argument("--platform", default="CUDA")
    parser.add_argument("--device", default="0")
    
    args = parser.parse_args()
    
    # Setup paths
    structures_dir = "/home/pjho3/projects/Drug/structures"
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Prepare receptor once
    cleaned_receptor = prepare_receptor(args.receptor, args.output_dir)
    
    # Define models to run
    models = []
    if args.model in ["a", "both"]:
        models.append(("model_a", f"model_a_tris_peg{args.peg}.sdf"))
    if args.model in ["b", "both"]:
        models.append(("model_b", f"model_b_arg_tris_peg{args.peg}.sdf"))
    
    # Run simulations
    results = {}
    for model_name, ligand_file in models:
        print("\n" + "=" * 60)
        print(f"Running {model_name.upper()} (PEG{args.peg})")
        print("=" * 60)
        
        ligand_sdf = os.path.join(structures_dir, ligand_file)
        if not os.path.exists(ligand_sdf):
            print(f"❌ Ligand file not found: {ligand_sdf}")
            continue
        
        # Create model output directory
        model_dir = os.path.join(args.output_dir, model_name)
        os.makedirs(model_dir, exist_ok=True)
        
        # Place ligand
        docked_sdf = os.path.join(model_dir, "docked_ligand.sdf")
        place_ligand_at_entrance(ligand_sdf, cleaned_receptor, docked_sdf)
        
        # Run MD
        output_prefix = os.path.join(model_dir, f"prod_{model_name}")
        traj, final, log = run_md_simulation(
            receptor_pdb=cleaned_receptor,
            ligand_sdf=docked_sdf,
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
