#!/usr/bin/env python
"""
Phase 4: Glycan Penetration Test
=================================
Compare drug accessibility between:
- Naked GLUT1 (cancer cell model)
- Glycosylated GLUT1 (normal cell model - RBC, Endothelial)

Hypothesis: Large tripod drug cannot penetrate glycan layer

Methods:
1. Place drug outside binding site (extracellular)
2. Run unbiased MD simulation
3. Measure COM distance to binding site over time
4. Analyze H-bonds with glycan vs protein
5. Calculate Solvent Accessible Surface Area (SAS)
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

import warnings
warnings.filterwarnings('ignore')

# Paths
STRUCTURES_DIR = "/home/pjho3/projects/Drug/structures/phase4"
LIGAND_SDF = "/home/pjho3/projects/Drug/structures/model_a_tris_peg2.sdf"
OUTPUT_DIR = "/home/pjho3/projects/Drug/results/phase4_glycan"

os.makedirs(OUTPUT_DIR, exist_ok=True)


def load_ligand(sdf_path):
    """Load ligand from SDF file"""
    print(f"\nLoading ligand: {sdf_path}")
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
    mol = next(suppl)
    if mol is None:
        raise ValueError(f"Could not load molecule from {sdf_path}")
    
    mol = Chem.AddHs(mol, addCoords=True)
    print(f"  Atoms: {mol.GetNumAtoms()}")
    return mol


def place_ligand_outside(receptor_pdb, ligand_mol, distance=30.0):
    """
    Place ligand outside the binding site, in extracellular space.
    The drug should approach from outside, potentially blocked by glycan.
    """
    print(f"\nPlacing ligand {distance}Å from binding site...")
    
    # Load receptor to find binding site
    pdb = PDBFile(receptor_pdb)
    positions = np.array([[p.x, p.y, p.z] for p in pdb.positions]) * 10  # nm to Å
    
    # Find protein center
    protein_center = np.mean(positions, axis=0)
    
    # Find residue 45 (glycosylation site) - this is near extracellular entrance
    res45_pos = None
    for atom in pdb.topology.atoms():
        if atom.residue.id == '45' and atom.name == 'CA':
            pos = pdb.positions[atom.index]
            res45_pos = np.array([pos.x, pos.y, pos.z]) * 10
            break
    
    if res45_pos is None:
        # Fallback: use approximate position
        res45_pos = protein_center + np.array([0, 30, 0])
    
    print(f"  Residue 45 position: ({res45_pos[0]:.1f}, {res45_pos[1]:.1f}, {res45_pos[2]:.1f})")
    
    # Direction from protein center to res45 (extracellular direction)
    direction = res45_pos - protein_center
    direction = direction / np.linalg.norm(direction)
    
    # Place ligand further out in this direction
    ligand_target = res45_pos + direction * distance
    
    print(f"  Ligand target: ({ligand_target[0]:.1f}, {ligand_target[1]:.1f}, {ligand_target[2]:.1f})")
    
    # Get ligand conformer and translate
    conf = ligand_mol.GetConformer()
    lig_positions = conf.GetPositions()
    lig_center = np.mean(lig_positions, axis=0)
    
    # Translate ligand to target position
    translation = ligand_target - lig_center
    for i in range(ligand_mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (pos.x + translation[0], 
                                  pos.y + translation[1], 
                                  pos.z + translation[2]))
    
    return ligand_mol, ligand_target


def run_penetration_simulation(receptor_pdb, ligand_mol, output_prefix,
                                model_name="naked",
                                total_ns=2.0, equil_ns=0.2,
                                temperature_k=310.0, dt_fs=2.0, save_ps=10.0):
    """
    Run MD simulation to test drug penetration.
    """
    print("\n" + "=" * 60)
    print(f"Running Penetration Test: {model_name.upper()} GLUT1")
    print("=" * 60)
    
    t0 = time.time()
    
    out_dir = os.path.join(OUTPUT_DIR, model_name)
    os.makedirs(out_dir, exist_ok=True)
    
    traj_dcd = os.path.join(out_dir, f"penetration_{model_name}.dcd")
    final_pdb = os.path.join(out_dir, f"penetration_{model_name}_final.pdb")
    log_csv = os.path.join(out_dir, f"penetration_{model_name}_log.csv")
    chk_path = os.path.join(out_dir, f"penetration_{model_name}.chk")
    
    # Prepare ligand for OpenMM
    print("   Preparing ligand...")
    ligand_off = Molecule.from_rdkit(ligand_mol, allow_undefined_stereo=True)
    ligand_off.assign_partial_charges(partial_charge_method="gasteiger")
    
    cache_file = os.path.join(out_dir, "gaff_cache.json")
    gaff = GAFFTemplateGenerator(molecules=[ligand_off], cache=cache_file)
    
    # Save ligand as PDB
    temp_lig_pdb = os.path.join(out_dir, "temp_lig.pdb")
    Chem.MolToPDBFile(ligand_mol, temp_lig_pdb)
    
    # Merge receptor and ligand
    print("   Merging complex...")
    with open(receptor_pdb, 'r') as f:
        lines_prot = [l for l in f.readlines() 
                      if not l.startswith("END") and not l.startswith("CONECT")]
    
    with open(temp_lig_pdb, 'r') as f:
        lines_lig = [l for l in f.readlines() 
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
    
    # Add solvent with larger box for drug movement
    print("   Adding solvent (larger box for penetration test)...")
    modeller.addSolvent(
        forcefield,
        padding=1.5 * nanometer,  # Larger padding
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
    print("\nComparing drug accessibility:")
    print("  - Naked GLUT1 (cancer cell)")
    print("  - Glycosylated GLUT1 (normal cell)")
    
    # Load ligand (TRIS model - the selected drug)
    ligand_mol = load_ligand(LIGAND_SDF)
    
    # Test 1: Naked GLUT1
    print("\n" + "=" * 70)
    print("TEST 1: NAKED GLUT1 (Cancer Cell Model)")
    print("=" * 70)
    
    naked_pdb = os.path.join(STRUCTURES_DIR, "glut1_naked.pdb")
    ligand_naked, _ = place_ligand_outside(naked_pdb, Chem.RWMol(ligand_mol), distance=25.0)
    
    naked_dcd, naked_final = run_penetration_simulation(
        naked_pdb, ligand_naked,
        output_prefix="naked",
        model_name="naked",
        total_ns=2.0,
        equil_ns=0.2
    )
    
    # Test 2: Glycosylated GLUT1
    print("\n" + "=" * 70)
    print("TEST 2: GLYCOSYLATED GLUT1 (Normal Cell Model)")
    print("=" * 70)
    
    glyco_pdb = os.path.join(STRUCTURES_DIR, "glut1_glycosylated.pdb")
    ligand_glyco, _ = place_ligand_outside(glyco_pdb, Chem.RWMol(ligand_mol), distance=25.0)
    
    glyco_dcd, glyco_final = run_penetration_simulation(
        glyco_pdb, ligand_glyco,
        output_prefix="glycosylated",
        model_name="glycosylated",
        total_ns=2.0,
        equil_ns=0.2
    )
    
    print("\n" + "=" * 70)
    print("Penetration Tests Complete!")
    print("=" * 70)
    print("\nNext: Run analysis script to compare results")


if __name__ == "__main__":
    main()
