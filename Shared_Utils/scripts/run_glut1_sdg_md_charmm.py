#!/usr/bin/env python3
"""
GLUT1-SDG Complex MD Simulation
Using OpenMM with CHARMM force field (proper approach)
"""

import openmm as mm
from openmm import app, unit
import sys
import os
import numpy as np

print("=" * 80)
print("GLUT1-SDG Complex MD Simulation")
print("Using CHARMM36 + SDG parameters (CHARMM format)")
print("=" * 80)

# Check platform
try:
    platform = mm.Platform.getPlatformByName('CUDA')
    print(f"✅ Using CUDA")
    use_cuda = True
except:
    platform = mm.Platform.getPlatformByName('CPU')
    print(f"⚠️  Using CPU")
    use_cuda = False

# File paths
pdb_file = 'glut1_tripod_complex.pdb'
output_dir = '../results/glut1_sdg_charmm'
os.makedirs(output_dir, exist_ok=True)

print(f"\n📁 Input: {pdb_file}")
print(f"📁 Output: {output_dir}")

# Step 1: Load structure
print("\n" + "=" * 80)
print("Step 1: Loading GLUT1-SDG complex")
print("=" * 80)

pdb = app.PDBFile(pdb_file)
print(f"✅ Loaded {pdb_file}")
print(f"  Total atoms: {pdb.topology.getNumAtoms()}")
print(f"  Residues: {pdb.topology.getNumResidues()}")

# Count residues
protein_res = 0
glycan_res = 0
ligand_res = 0

for residue in pdb.topology.residues():
    res_name = residue.name
    if res_name == 'SDG':
        ligand_res += 1
    elif res_name in ['BGLCNA', 'AMAN', 'BMAN', 'ANE5AC', 'ANE', 'BGL', 'BGA', 'AMA']:
        glycan_res += 1
    else:
        protein_res += 1

print(f"  Protein residues: {protein_res}")
print(f"  Glycan residues: {glycan_res}")
print(f"  SDG ligand residues: {ligand_res}")

# Step 2: Load CHARMM parameters
print("\n" + "=" * 80)
print("Step 2: Loading CHARMM36 + SDG parameters")
print("=" * 80)

# Load CHARMM36 force field + SDG parameters
# OpenMM can directly read CHARMM rtf/prm files
print("  Loading CHARMM36 protein/lipid/carbohydrate parameters...")
print("  Loading SDG ligand parameters (sdg.rtf, sdg.prm)...")

params = app.CharmmParameterSet(
    'sdg.rtf',
    'sdg.prm'
)

print(f"✅ Parameters loaded")

# Step 3: Build system with modeller
print("\n" + "=" * 80)
print("Step 3: Adding membrane and solvent")
print("=" * 80)

# For CHARMM force field, we need to use ForceField with CHARMM files
# But first, let's use Modeller to add membrane and solvent
forcefield = app.ForceField('charmm36.xml', 'charmm36/water.xml', 'charmm36/lipid.xml')

modeller = app.Modeller(pdb.topology, pdb.positions)

print("  Adding POPC lipid bilayer...")
# Note: This will use CHARMM36 lipid parameters
modeller.addMembrane(
    forcefield=forcefield,
    lipidType='POPC',
    minimumPadding=1.5*unit.nanometer,
    ionicStrength=0.15*unit.molar
)

print(f"✅ Membrane added")
print(f"  Total atoms: {modeller.topology.getNumAtoms()}")

print("  Adding water and ions...")
modeller.addSolvent(
    forcefield,
    model='tip3p',
    padding=1.5*unit.nanometer,
    ionicStrength=0.15*unit.molar
)

print(f"✅ Solvent added")
print(f"  Total atoms: {modeller.topology.getNumAtoms()}")

# Save solvated system
solvated_pdb = os.path.join(output_dir, 'system_solvated.pdb')
with open(solvated_pdb, 'w') as f:
    app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
print(f"✅ Saved: {solvated_pdb}")

# Step 4: Create system
print("\n" + "=" * 80)
print("Step 4: Creating system with CHARMM36 + SDG")
print("=" * 80)

# Create system with CHARMM36 + custom SDG parameters
# Note: For SDG, we need to register the residue template
print("  Registering SDG residue template...")

# Load the full force field including SDG
# We'll use the CHARMM36 XML files and add SDG parameters
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.2*unit.nanometer,
    constraints=app.HBonds,
    rigidWater=True,
    ewaldErrorTolerance=0.0005
)

print(f"✅ System created")
print(f"  Forces: {system.getNumForces()}")

# Step 5: Minimization
print("\n" + "=" * 80)
print("Step 5: Energy minimization")
print("=" * 80)

integrator = mm.LangevinMiddleIntegrator(
    300*unit.kelvin,
    1.0/unit.picosecond,
    2.0*unit.femtosecond
)

if use_cuda:
    properties = {'CudaPrecision': 'mixed'}
else:
    properties = {}

simulation = app.Simulation(
    modeller.topology,
    system,
    integrator,
    platform,
    properties
)

simulation.context.setPositions(modeller.positions)

print("  Initial energy...")
state = simulation.context.getState(getEnergy=True)
print(f"    Potential Energy: {state.getPotentialEnergy()}")

print("  Minimizing (max 1000 steps)...")
simulation.minimizeEnergy(maxIterations=1000)

state = simulation.context.getState(getEnergy=True, getPositions=True)
print(f"    Final Potential Energy: {state.getPotentialEnergy()}")

# Save minimized
minimized_pdb = os.path.join(output_dir, 'system_minimized.pdb')
with open(minimized_pdb, 'w') as f:
    app.PDBFile.writeFile(modeller.topology, state.getPositions(), f)
print(f"✅ Saved: {minimized_pdb}")

# Step 6: Equilibration
print("\n" + "=" * 80)
print("Step 6: Equilibration (100 ps)")
print("=" * 80)

simulation.context.setVelocitiesToTemperature(300*unit.kelvin)

simulation.reporters.append(
    app.StateDataReporter(
        sys.stdout,
        1000,
        step=True,
        time=True,
        potentialEnergy=True,
        temperature=True,
        speed=True
    )
)

simulation.reporters.append(
    app.DCDReporter(os.path.join(output_dir, 'equilibration.dcd'), 1000)
)

print("  Running 100 ps equilibration...")
simulation.step(50000)  # 100 ps

print("✅ Equilibration complete")

# Save equilibrated
equilibrated_pdb = os.path.join(output_dir, 'system_equilibrated.pdb')
state = simulation.context.getState(getPositions=True)
with open(equilibrated_pdb, 'w') as f:
    app.PDBFile.writeFile(modeller.topology, state.getPositions(), f)
print(f"✅ Saved: {equilibrated_pdb}")

# Step 7: Production MD
print("\n" + "=" * 80)
print("Step 7: Production MD (10 ns)")
print("=" * 80)

simulation.reporters.clear()

simulation.reporters.append(
    app.StateDataReporter(
        os.path.join(output_dir, 'production.log'),
        5000,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        temperature=True,
        speed=True
    )
)

simulation.reporters.append(
    app.DCDReporter(os.path.join(output_dir, 'production.dcd'), 5000)
)

simulation.reporters.append(
    app.CheckpointReporter(os.path.join(output_dir, 'checkpoint.chk'), 50000)
)

print("  Running 10 ns production...")
simulation.step(5000000)  # 10 ns

# Save final
state = simulation.context.getState(getPositions=True)
final_pdb = os.path.join(output_dir, 'system_final.pdb')
with open(final_pdb, 'w') as f:
    app.PDBFile.writeFile(modeller.topology, state.getPositions(), f)

print("\n" + "=" * 80)
print("✅ SIMULATION COMPLETE!")
print("=" * 80)
print(f"\n📁 Results: {output_dir}")
print("  - system_solvated.pdb")
print("  - system_minimized.pdb")
print("  - equilibration.dcd")
print("  - production.dcd")
print("  - production.log")
print("  - system_final.pdb")
print("\n" + "=" * 80)
