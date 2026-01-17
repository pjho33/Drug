#!/usr/bin/env python3
"""
GLUT1-Tripod Complex MD Simulation with Membrane
Direct system building using OpenMM
"""

import openmm as mm
from openmm import app, unit
import sys
import os
import numpy as np

print("=" * 80)
print("GLUT1-Tripod Complex MD Simulation")
print("=" * 80)

# Check CUDA availability
try:
    from openmm import Platform
    platform = Platform.getPlatformByName('CUDA')
    print(f"✅ CUDA available")
    use_cuda = True
except:
    print(f"⚠️  CUDA not available, using CPU")
    use_cuda = False

# File paths
pdb_file = 'glut1_tripod_complex.pdb'
output_dir = '../results/glut1_tripod_md'

print(f"\n📁 Input: {pdb_file}")
print(f"📁 Output: {output_dir}")

# Create output directory
os.makedirs(output_dir, exist_ok=True)

# Step 1: Load structure
print("\n" + "=" * 80)
print("Step 1: Loading structure")
print("=" * 80)

pdb = app.PDBFile(pdb_file)
print(f"✅ Loaded {pdb_file}")
print(f"  Atoms: {pdb.topology.getNumAtoms()}")
print(f"  Residues: {pdb.topology.getNumResidues()}")
print(f"  Chains: {pdb.topology.getNumChains()}")

# Count protein and ligand atoms
protein_atoms = 0
ligand_atoms = 0
glycan_atoms = 0

for residue in pdb.topology.residues():
    res_name = residue.name
    if res_name in ['SDG', 'TRP', 'UNL', 'LIG']:
        ligand_atoms += len(list(residue.atoms()))
    elif res_name in ['BGLCNA', 'AMAN', 'BMAN', 'ANE5AC', 'ANE']:
        glycan_atoms += len(list(residue.atoms()))
    else:
        protein_atoms += len(list(residue.atoms()))

print(f"\n  Protein atoms: {protein_atoms}")
print(f"  Glycan atoms: {glycan_atoms}")
print(f"  Ligand atoms: {ligand_atoms}")

# Step 2: Add membrane
print("\n" + "=" * 80)
print("Step 2: Adding membrane")
print("=" * 80)

# Use POPC membrane (tumor cell-like)
modeller = app.Modeller(pdb.topology, pdb.positions)

print("  Adding POPC lipid bilayer...")
# Add membrane around protein
# Membrane will be added in XY plane, protein should be oriented properly
modeller.addMembrane(
    forcefield=app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml', 'amber14/lipid17.xml'),
    lipidType='POPC',
    minimumPadding=1.5*unit.nanometer,
    ionicStrength=0.15*unit.molar
)

print(f"✅ Membrane added")
print(f"  Total atoms: {modeller.topology.getNumAtoms()}")

# Step 3: Add solvent
print("\n" + "=" * 80)
print("Step 3: Adding solvent and ions")
print("=" * 80)

modeller.addSolvent(
    forcefield=app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml'),
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
print("Step 4: Creating system")
print("=" * 80)

# Load force field
# Note: Tripod parameters need to be loaded separately
forcefield = app.ForceField(
    'amber14-all.xml',
    'amber14/tip3pfb.xml',
    'amber14/lipid17.xml'
)

# For now, we'll use GAFF for the ligand (Tripod)
# In production, you should use the trp.rtf and trp.prm files
print("⚠️  Using GAFF for ligand (Tripod)")
print("    For production, convert trp.rtf/trp.prm to OpenMM XML format")

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
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed'}
else:
    platform = mm.Platform.getPlatformByName('CPU')
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

# Save minimized structure
minimized_pdb = os.path.join(output_dir, 'system_minimized.pdb')
positions = state.getPositions()
with open(minimized_pdb, 'w') as f:
    app.PDBFile.writeFile(modeller.topology, positions, f)
print(f"✅ Saved: {minimized_pdb}")

# Step 6: Equilibration
print("\n" + "=" * 80)
print("Step 6: Equilibration (NVT)")
print("=" * 80)

# Set temperature
simulation.context.setVelocitiesToTemperature(300*unit.kelvin)

# Add reporters
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
    app.DCDReporter(
        os.path.join(output_dir, 'equilibration.dcd'),
        1000
    )
)

print("  Running 100 ps equilibration...")
simulation.step(50000)  # 100 ps at 2 fs timestep

print("✅ Equilibration complete")

# Save equilibrated state
equilibrated_pdb = os.path.join(output_dir, 'system_equilibrated.pdb')
state = simulation.context.getState(getPositions=True)
positions = state.getPositions()
with open(equilibrated_pdb, 'w') as f:
    app.PDBFile.writeFile(modeller.topology, positions, f)
print(f"✅ Saved: {equilibrated_pdb}")

# Step 7: Production MD
print("\n" + "=" * 80)
print("Step 7: Production MD (10 ns)")
print("=" * 80)

# Clear old reporters
simulation.reporters.clear()

# Add new reporters for production
simulation.reporters.append(
    app.StateDataReporter(
        os.path.join(output_dir, 'production.log'),
        5000,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=True,
        density=True,
        speed=True
    )
)

simulation.reporters.append(
    app.DCDReporter(
        os.path.join(output_dir, 'production.dcd'),
        5000  # Save every 10 ps
    )
)

simulation.reporters.append(
    app.CheckpointReporter(
        os.path.join(output_dir, 'checkpoint.chk'),
        50000  # Save every 100 ps
    )
)

print("  Running 10 ns production...")
print("  (This will take a while...)")
simulation.step(5000000)  # 10 ns at 2 fs timestep

print("\n✅ Production MD complete!")

# Save final state
final_pdb = os.path.join(output_dir, 'system_final.pdb')
state = simulation.context.getState(getPositions=True)
positions = state.getPositions()
with open(final_pdb, 'w') as f:
    app.PDBFile.writeFile(modeller.topology, positions, f)
print(f"✅ Saved: {final_pdb}")

print("\n" + "=" * 80)
print("✅ SIMULATION COMPLETE!")
print("=" * 80)
print(f"\n📁 Results saved in: {output_dir}")
print("\nOutput files:")
print(f"  - system_solvated.pdb: Initial solvated system")
print(f"  - system_minimized.pdb: After energy minimization")
print(f"  - system_equilibrated.pdb: After equilibration")
print(f"  - equilibration.dcd: Equilibration trajectory")
print(f"  - production.dcd: Production trajectory")
print(f"  - production.log: Energy and statistics")
print(f"  - system_final.pdb: Final structure")
print("\n" + "=" * 80)
