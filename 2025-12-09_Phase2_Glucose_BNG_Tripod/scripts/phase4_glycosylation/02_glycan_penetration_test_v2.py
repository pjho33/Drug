#!/usr/bin/env python
"""
Phase 4: Glycan Penetration Test (CHARMM36 + SystemGenerator)
==============================================================
Compare drug accessibility between:
- Naked GLUT1 (cancer cell model)
- Glycosylated GLUT1 (normal cell model - RBC, Endothelial)

Uses SystemGenerator to handle glycan residues (NAG, MAN, BMA)
with CHARMM36 forcefield + GAFF for ligand.
"""

import os
import sys
import time
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem

from openmm import Platform, MonteCarloBarostat, LangevinMiddleIntegrator
from openmm.app import (
    DCDReporter, PDBFile, Simulation,
    StateDataReporter, CheckpointReporter, PME, HBonds, Modeller
)
from openmm.unit import bar, kelvin, molar, nanometer, picosecond, picoseconds

from openmmforcefields.generators import SystemGenerator
from openff.toolkit.topology import Molecule

from pdbfixer import PDBFixer

import warnings
warnings.filterwarnings('ignore')

# Paths
STRUCTURES_DIR = "/home/pjho3/projects/Drug/structures/phase4"
LIGAND_SDF = "/home/pjho3/projects/Drug/structures/model_a_tris_peg2.sdf"
OUTPUT_DIR = "/home/pjho3/projects/Drug/results/phase4_glycan"

os.makedirs(OUTPUT_DIR, exist_ok=True)


def load_and_prepare_ligand(sdf_path):
    """Load ligand and prepare for OpenMM"""
    print(f"\n📦 Loading ligand: {sdf_path}")
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
    mol = next(suppl)
    if mol is None:
        raise ValueError(f"Could not load molecule from {sdf_path}")
    
    mol = Chem.AddHs(mol, addCoords=True)
    print(f"   Atoms: {mol.GetNumAtoms()}")
    
    # Create OpenFF molecule
    ligand_off = Molecule.from_rdkit(mol, allow_undefined_stereo=True)
    
    return mol, ligand_off


def place_ligand_outside(receptor_pdb_path, ligand_mol, distance=30.0):
    """Place ligand outside binding site in extracellular space"""
    print(f"\n📍 Placing ligand {distance}Å from protein...")
    
    # Load receptor
    pdb = PDBFile(receptor_pdb_path)
    positions = np.array([[p.x, p.y, p.z] for p in pdb.positions]) * 10  # nm to Å
    
    # Find protein center and residue 45 (glycosylation site)
    protein_center = np.mean(positions, axis=0)
    
    res45_pos = None
    for atom in pdb.topology.atoms():
        if atom.residue.id == '45' and atom.name == 'CA':
            pos = pdb.positions[atom.index]
            res45_pos = np.array([pos.x, pos.y, pos.z]) * 10
            break
    
    if res45_pos is None:
        res45_pos = protein_center + np.array([0, 30, 0])
    
    # Direction outward from protein
    direction = res45_pos - protein_center
    direction = direction / np.linalg.norm(direction)
    
    # Target position for ligand
    ligand_target = res45_pos + direction * distance
    print(f"   Target position: ({ligand_target[0]:.1f}, {ligand_target[1]:.1f}, {ligand_target[2]:.1f})")
    
    # Translate ligand
    conf = ligand_mol.GetConformer()
    lig_positions = conf.GetPositions()
    lig_center = np.mean(lig_positions, axis=0)
    translation = ligand_target - lig_center
    
    for i in range(ligand_mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (pos.x + translation[0],
                                  pos.y + translation[1],
                                  pos.z + translation[2]))
    
    return ligand_mol


def prepare_receptor(pdb_path, output_dir):
    """Prepare receptor with PDBFixer"""
    print(f"\n🔧 Preparing receptor: {pdb_path}")
    
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    
    prepared_pdb = os.path.join(output_dir, "receptor_prepared.pdb")
    with open(prepared_pdb, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    
    print(f"   Prepared receptor saved: {prepared_pdb}")
    return prepared_pdb


def run_penetration_simulation(receptor_pdb, ligand_mol, ligand_off, output_dir,
                                model_name="naked",
                                total_ns=2.0, equil_ns=0.2,
                                temperature_k=310.0, dt_fs=2.0, save_ps=10.0):
    """Run MD simulation with SystemGenerator for glycan support"""
    print("\n" + "=" * 60)
    print(f"🚀 Running Penetration Test: {model_name.upper()}")
    print("=" * 60)
    
    t0 = time.time()
    
    out_dir = os.path.join(output_dir, model_name)
    os.makedirs(out_dir, exist_ok=True)
    
    traj_dcd = os.path.join(out_dir, f"penetration_{model_name}.dcd")
    final_pdb = os.path.join(out_dir, f"penetration_{model_name}_final.pdb")
    log_csv = os.path.join(out_dir, f"penetration_{model_name}_log.csv")
    chk_path = os.path.join(out_dir, f"penetration_{model_name}.chk")
    
    # Save ligand as PDB
    temp_lig_pdb = os.path.join(out_dir, "temp_lig.pdb")
    Chem.MolToPDBFile(ligand_mol, temp_lig_pdb)
    
    # Merge receptor and ligand PDBs
    print("   Merging receptor and ligand...")
    with open(receptor_pdb, 'r') as f:
        receptor_lines = [l for l in f.readlines() 
                         if l.startswith("ATOM") or l.startswith("HETATM")]
    
    with open(temp_lig_pdb, 'r') as f:
        ligand_lines = [l for l in f.readlines()
                       if l.startswith("ATOM") or l.startswith("HETATM")]
    
    merged_pdb = os.path.join(out_dir, "complex_merged.pdb")
    with open(merged_pdb, 'w') as f:
        f.writelines(receptor_lines)
        f.write("TER\n")
        f.writelines(ligand_lines)
        f.write("END\n")
    
    # Load merged complex
    print("   Loading merged complex...")
    pdb = PDBFile(merged_pdb)
    modeller = Modeller(pdb.topology, pdb.positions)
    
    # Create SystemGenerator with CHARMM36 + GAFF for ligand
    print("   Creating SystemGenerator (CHARMM36 + GAFF)...")
    system_generator = SystemGenerator(
        forcefields=['charmm36.xml', 'charmm36/water.xml'],
        small_molecule_forcefield='gaff-2.11',
        molecules=[ligand_off],
        cache=os.path.join(out_dir, 'gaff_cache.json')
    )
    
    # Add solvent
    print("   Adding solvent...")
    modeller.addSolvent(
        system_generator.forcefield,
        padding=1.2 * nanometer,
        ionicStrength=0.15 * molar
    )
    
    print(f"   Total atoms: {modeller.topology.getNumAtoms()}")
    
    # Create system
    print("   Creating system...")
    system = system_generator.create_system(modeller.topology)
    
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
    
    # Minimize
    print("   Minimizing energy (phase 1)...")
    simulation.minimizeEnergy(maxIterations=5000, tolerance=100)
    print("   Minimizing energy (phase 2)...")
    simulation.minimizeEnergy(maxIterations=10000, tolerance=10)
    
    # Gradual heating
    print("   Gradual heating...")
    simulation.context.setVelocitiesToTemperature(50 * kelvin)
    for temp in [50, 100, 150, 200, 250, 300, temperature_k]:
        integrator.setTemperature(temp * kelvin)
        simulation.step(1000)
        print(f"     Heated to {temp} K")
    
    # Calculate steps
    total_steps = int(total_ns * 1e6 / dt_fs)
    equil_steps = int(equil_ns * 1e6 / dt_fs)
    save_freq = int(save_ps * 1000 / dt_fs)
    prod_steps = total_steps - equil_steps
    
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
            print(f"     Equilibration: {(i + steps_to_run) / equil_steps * 100:.0f}%")
    
    # Production reporters
    simulation.reporters.append(DCDReporter(traj_dcd, save_freq))
    simulation.reporters.append(StateDataReporter(
        log_csv, save_freq,
        step=True, potentialEnergy=True, temperature=True,
        volume=True, speed=True
    ))
    simulation.reporters.append(CheckpointReporter(chk_path, 5000))
    
    # Production
    print(f"\n   Running production for {total_ns - equil_ns} ns...")
    chunk_size = 50000
    steps_done = 0
    while steps_done < prod_steps:
        steps_to_run = min(chunk_size, prod_steps - steps_done)
        simulation.step(steps_to_run)
        steps_done += steps_to_run
        print(f"     Production: {steps_done / prod_steps * 100:.1f}%")
    
    # Save final structure
    state = simulation.context.getState(getPositions=True)
    with open(final_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, state.getPositions(), f)
    
    elapsed = time.time() - t0
    print(f"\n✅ Simulation complete!")
    print(f"   Time: {elapsed:.1f} seconds")
    print(f"   Trajectory: {traj_dcd}")
    
    return traj_dcd, final_pdb


def main():
    print("=" * 70)
    print("Phase 4: Glycan Penetration Test")
    print("=" * 70)
    print("\nHypothesis: Large drug cannot penetrate glycan layer on normal cells")
    print("\nModels:")
    print("  - Naked GLUT1 (cancer cell) - no glycan barrier")
    print("  - Glycosylated GLUT1 (normal cell) - Man5 glycan at Asn45")
    
    # Load ligand
    ligand_mol, ligand_off = load_and_prepare_ligand(LIGAND_SDF)
    
    # Test 1: Naked GLUT1
    print("\n" + "=" * 70)
    print("TEST 1: NAKED GLUT1 (Cancer Cell Model)")
    print("=" * 70)
    
    naked_pdb = os.path.join(STRUCTURES_DIR, "glut1_naked_asn45.pdb")
    naked_out_dir = os.path.join(OUTPUT_DIR, "naked")
    os.makedirs(naked_out_dir, exist_ok=True)
    
    # Prepare receptor
    naked_prepared = prepare_receptor(naked_pdb, naked_out_dir)
    
    # Place ligand outside
    ligand_naked = place_ligand_outside(naked_prepared, Chem.RWMol(ligand_mol), distance=25.0)
    
    # Run simulation
    naked_dcd, naked_final = run_penetration_simulation(
        naked_prepared, ligand_naked, ligand_off,
        OUTPUT_DIR, model_name="naked",
        total_ns=2.0, equil_ns=0.2
    )
    
    # Test 2: Glycosylated GLUT1
    print("\n" + "=" * 70)
    print("TEST 2: GLYCOSYLATED GLUT1 (Normal Cell Model)")
    print("=" * 70)
    
    glyco_pdb = os.path.join(STRUCTURES_DIR, "glut1_glycosylated_man5.pdb")
    glyco_out_dir = os.path.join(OUTPUT_DIR, "glycosylated")
    os.makedirs(glyco_out_dir, exist_ok=True)
    
    # Prepare receptor (with glycan)
    glyco_prepared = prepare_receptor(glyco_pdb, glyco_out_dir)
    
    # Place ligand outside
    ligand_glyco = place_ligand_outside(glyco_prepared, Chem.RWMol(ligand_mol), distance=25.0)
    
    # Run simulation
    glyco_dcd, glyco_final = run_penetration_simulation(
        glyco_prepared, ligand_glyco, ligand_off,
        OUTPUT_DIR, model_name="glycosylated",
        total_ns=2.0, equil_ns=0.2
    )
    
    print("\n" + "=" * 70)
    print("✅ All Penetration Tests Complete!")
    print("=" * 70)
    print("\nRun analysis script to compare results:")
    print("  python scripts/phase4_glycosylation/03_analyze_penetration.py")


if __name__ == "__main__":
    main()
