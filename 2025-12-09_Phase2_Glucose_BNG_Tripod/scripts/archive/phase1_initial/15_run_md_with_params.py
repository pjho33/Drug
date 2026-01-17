#!/usr/bin/env python
"""
Phase 2: MD Simulation with Pre-parameterized Ligands
=====================================================
Uses mol2/frcmod files generated from Gasteiger charges
to avoid AM1BCC timeout issues with large molecules.
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
from Bio.PDB import PDBParser


def prepare_receptor(pdb_path, output_dir):
    """Prepare GLUT1 receptor using PDBFixer"""
    print("=" * 60)
    print("Step 1: Preparing Receptor")
    print("=" * 60)
    
    cleaned_path = os.path.join(output_dir, "cleaned_receptor.pdb")
    
    # Check if already prepared
    if os.path.exists(cleaned_path):
        print(f"   Using existing: {cleaned_path}")
        return cleaned_path
    
    fixer = PDBFixer(filename=pdb_path)
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    
    with open(cleaned_path, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    
    print(f"   ✅ Saved: {cleaned_path}")
    return cleaned_path


def get_binding_site_center(receptor_pdb):
    """Get GLUT1 binding site center"""
    BINDING_RESIDUES = [34, 161, 282, 283, 288, 317, 388, 412]
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('receptor', receptor_pdb)
    
    all_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    all_coords.append(residue['CA'].get_coord())
    
    all_coords = np.array(all_coords)
    center = all_coords.mean(axis=0)
    z_max = all_coords[:, 2].max()
    
    return center, z_max


def place_ligand(sdf_file, receptor_pdb, output_sdf):
    """Place ligand at GLUT1 channel entrance"""
    print("\n" + "=" * 60)
    print("Step 2: Placing Ligand")
    print("=" * 60)
    
    # Load from SDF (more reliable than mol2 with RDKit)
    suppl = Chem.SDMolSupplier(sdf_file, removeHs=False)
    mol = next(suppl)
    if mol is None:
        print(f"   ❌ Failed to load {sdf_file}")
        return None
    
    conf = mol.GetConformer()
    center, z_max = get_binding_site_center(receptor_pdb)
    
    print(f"   Receptor center: ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})")
    print(f"   Z_max: {z_max:.1f} Å")
    
    # Center ligand
    lig_coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    lig_center = lig_coords.mean(axis=0)
    
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Geometry.Point3D(
            float(pos.x - lig_center[0]),
            float(pos.y - lig_center[1]),
            float(pos.z - lig_center[2])))
    
    # Check arm orientation and flip if needed
    lig_coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    distances = np.linalg.norm(lig_coords, axis=1)
    arm_z_avg = lig_coords[distances > np.percentile(distances, 90), 2].mean()
    
    if arm_z_avg > 0:
        print("   🔄 Flipping ligand")
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, Geometry.Point3D(float(pos.x), float(-pos.y), float(-pos.z)))
    
    # Position above channel
    target_z = z_max + 15
    target_pos = np.array([center[0], center[1], target_z])
    
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Geometry.Point3D(
            float(pos.x + target_pos[0]),
            float(pos.y + target_pos[1]),
            float(pos.z + target_pos[2])))
    
    final_coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    print(f"   Ligand Z: {final_coords[:, 2].mean():.1f} Å")
    
    # Save as SDF
    writer = Chem.SDWriter(output_sdf)
    writer.write(mol)
    writer.close()
    print(f"   ✅ Saved: {output_sdf}")
    
    return output_sdf


def run_md(receptor_pdb, ligand_sdf, output_prefix,
           total_ns=100.0, equil_ns=1.0,
           temperature_k=310.0, dt_fs=2.0, save_ps=10.0,
           platform_name="CUDA"):
    """Run MD simulation"""
    print("\n" + "=" * 60)
    print("Step 3: MD Simulation")
    print("=" * 60)
    
    t0 = time.time()
    out_dir = os.path.dirname(output_prefix)
    os.makedirs(out_dir, exist_ok=True)
    
    traj_dcd = output_prefix + ".dcd"
    final_pdb = output_prefix + "_final.pdb"
    log_csv = output_prefix + "_log.csv"
    chk_path = output_prefix + ".chk"
    
    # Load ligand
    print("   Loading ligand...")
    suppl = Chem.SDMolSupplier(ligand_sdf, removeHs=False)
    ligand_rdkit = next(suppl)
    if ligand_rdkit is None:
        raise ValueError("Failed to load ligand")
    
    # Create OpenFF molecule with Gasteiger charges pre-assigned
    ligand_off = Molecule.from_rdkit(ligand_rdkit, allow_undefined_stereo=True)
    
    # Assign Gasteiger charges (this is fast)
    print("   Assigning Gasteiger charges...")
    ligand_off.assign_partial_charges(partial_charge_method="gasteiger")
    
    # Setup GAFF with pre-charged molecule
    cache_file = os.path.join(out_dir, "gaff_cache.json")
    gaff = GAFFTemplateGenerator(molecules=[ligand_off], cache=cache_file)
    
    # Merge complex
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
    
    pdb = PDBFile(merged_pdb)
    modeller = Modeller(pdb.topology, pdb.positions)
    
    # Add solvent
    print("   Adding solvent (this may take a few minutes)...")
    modeller.addSolvent(
        forcefield,
        padding=2.0 * nanometer,
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
    system.addForce(MonteCarloBarostat(1.0 * bar, temperature_k * kelvin))
    
    integrator = LangevinMiddleIntegrator(
        temperature_k * kelvin, 1 / picosecond, (dt_fs / 1000.0) * picoseconds
    )
    
    # Platform
    try:
        platform = Platform.getPlatformByName(platform_name)
        prop = {"DeviceIndex": "0", "Precision": "mixed"}
        print(f"   Platform: {platform_name}")
    except:
        platform = Platform.getPlatformByName("CPU")
        prop = {}
        print("   Platform: CPU (fallback)")
    
    simulation = Simulation(modeller.topology, system, integrator, platform, prop)
    simulation.context.setPositions(modeller.positions)
    
    # Minimize
    print("   Minimizing energy...")
    simulation.minimizeEnergy()
    
    simulation.context.setVelocitiesToTemperature(temperature_k * kelvin)
    
    # Steps
    eq_steps = int((equil_ns * 1_000_000.0) / dt_fs)
    prod_steps = int((total_ns * 1_000_000.0) / dt_fs)
    report_steps = max(1, int((save_ps * 1000.0) / dt_fs))
    
    print(f"   Equilibration: {equil_ns} ns")
    print(f"   Production: {total_ns} ns")
    
    # Reporters
    simulation.reporters.append(DCDReporter(traj_dcd, report_steps))
    simulation.reporters.append(StateDataReporter(
        log_csv, report_steps, step=True, potentialEnergy=True,
        temperature=True, volume=True, speed=True
    ))
    simulation.reporters.append(CheckpointReporter(chk_path, report_steps * 10))
    
    # Check for checkpoint
    if os.path.exists(chk_path):
        print(f"   📂 Loading checkpoint")
        simulation.loadCheckpoint(chk_path)
    else:
        if eq_steps > 0:
            print(f"   ⏳ Equilibrating...")
            simulation.step(eq_steps)
    
    # Production
    print(f"   ⏳ Running production...")
    simulation.step(prod_steps)
    
    # Save final
    state = simulation.context.getState(getPositions=True)
    with open(final_pdb, "w") as f:
        PDBFile.writeFile(simulation.topology, state.getPositions(), f)
    
    # Cleanup
    for f in [temp_lig_pdb, merged_pdb]:
        if os.path.exists(f):
            try:
                os.remove(f)
            except:
                pass
    
    dt = time.time() - t0
    print(f"\n✅ Complete! ({dt/3600:.1f} hours)")
    print(f"   Trajectory: {traj_dcd}")
    print(f"   Final: {final_pdb}")
    
    return traj_dcd, final_pdb, log_csv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=["a", "b", "both"], default="both")
    parser.add_argument("--total_ns", type=float, default=100.0)
    parser.add_argument("--equil_ns", type=float, default=1.0)
    parser.add_argument("--receptor", default="/home/pjho3/projects/Drug/raw_data/4PYP.pdb")
    parser.add_argument("--output_dir", default="/home/pjho3/projects/Drug/results/phase2_short_tripod")
    parser.add_argument("--platform", default="CUDA")
    
    args = parser.parse_args()
    
    params_dir = "/home/pjho3/projects/Drug/structures/parameterized"
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Prepare receptor
    cleaned_receptor = prepare_receptor(args.receptor, args.output_dir)
    
    # Models - use SDF files instead of mol2 (RDKit compatibility)
    structures_dir = "/home/pjho3/projects/Drug/structures"
    models = []
    if args.model in ["a", "both"]:
        models.append(("model_a", "model_a_tris_peg2.sdf"))
    if args.model in ["b", "both"]:
        models.append(("model_b", "model_b_arg_tris_peg2.sdf"))
    
    results = {}
    for model_name, sdf_file in models:
        print("\n" + "="*60)
        print(f"Running {model_name.upper()}")
        print("="*60)
        
        sdf_path = os.path.join(structures_dir, sdf_file)
        if not os.path.exists(sdf_path):
            print(f"❌ Not found: {sdf_path}")
            continue
        
        model_dir = os.path.join(args.output_dir, model_name)
        os.makedirs(model_dir, exist_ok=True)
        
        # Place ligand
        docked_sdf = os.path.join(model_dir, "docked_ligand.sdf")
        place_ligand(sdf_path, cleaned_receptor, docked_sdf)
        
        # Run MD
        output_prefix = os.path.join(model_dir, f"prod_{model_name}")
        traj, final, log = run_md(
            receptor_pdb=cleaned_receptor,
            ligand_sdf=docked_sdf,
            output_prefix=output_prefix,
            total_ns=args.total_ns,
            equil_ns=args.equil_ns,
            platform_name=args.platform
        )
        
        results[model_name] = {"trajectory": traj, "final_pdb": final, "log": log}
    
    # Save summary
    summary_path = os.path.join(args.output_dir, "simulation_summary.json")
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2)
    
    print("\n" + "="*60)
    print("Phase 2 Complete!")
    print("="*60)
    print(f"Results: {args.output_dir}")


if __name__ == "__main__":
    main()
