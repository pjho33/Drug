#!/usr/bin/env python3
"""
Simplified GLUT1-Tripod MD Simulation
Using OpenMM with automatic parameterization
"""

import openmm as mm
from openmm import app, unit
import sys
import os

print("=" * 80)
print("GLUT1-Tripod Complex MD Simulation (Simplified)")
print("=" * 80)

# Check platform
try:
    platform = mm.Platform.getPlatformByName('CUDA')
    print(f"✅ Using CUDA")
    use_cuda = True
except:
    platform = mm.Platform.getPlatformByName('CPU')
    print(f"⚠️  Using CPU (slower)")
    use_cuda = False

# Files
pdb_file = 'glut1_tripod_complex.pdb'
output_dir = '../results/glut1_tripod_simple'
os.makedirs(output_dir, exist_ok=True)

print(f"\n📁 Input: {pdb_file}")
print(f"📁 Output: {output_dir}")

# Load structure
print("\n" + "=" * 80)
print("Step 1: Loading structure")
print("=" * 80)

pdb = app.PDBFile(pdb_file)
print(f"✅ Loaded: {pdb.topology.getNumAtoms()} atoms")

# Create modeller
modeller = app.Modeller(pdb.topology, pdb.positions)

# Add solvent box (no membrane for now - simpler)
print("\n" + "=" * 80)
print("Step 2: Adding solvent box")
print("=" * 80)

forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

print("  Adding water box (1.5 nm padding)...")
modeller.addSolvent(
    forcefield,
    model='tip3p',
    padding=1.5*unit.nanometer,
    ionicStrength=0.15*unit.molar
)

print(f"✅ Total atoms: {modeller.topology.getNumAtoms()}")

# Save solvated
solvated_pdb = os.path.join(output_dir, 'solvated.pdb')
with open(solvated_pdb, 'w') as f:
    app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
print(f"✅ Saved: {solvated_pdb}")

# Create system
print("\n" + "=" * 80)
print("Step 3: Creating system")
print("=" * 80)

print("⚠️  Note: Tripod will use generic AMBER parameters")
print("    For production, use proper CHARMM parameters")

system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0*unit.nanometer,
    constraints=app.HBonds
)

print(f"✅ System created")

# Integrator
integrator = mm.LangevinMiddleIntegrator(
    300*unit.kelvin,
    1.0/unit.picosecond,
    2.0*unit.femtosecond
)

# Simulation
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

# Minimize
print("\n" + "=" * 80)
print("Step 4: Energy minimization")
print("=" * 80)

state = simulation.context.getState(getEnergy=True)
print(f"  Initial: {state.getPotentialEnergy()}")

simulation.minimizeEnergy(maxIterations=1000)

state = simulation.context.getState(getEnergy=True, getPositions=True)
print(f"  Final: {state.getPotentialEnergy()}")

minimized_pdb = os.path.join(output_dir, 'minimized.pdb')
with open(minimized_pdb, 'w') as f:
    app.PDBFile.writeFile(modeller.topology, state.getPositions(), f)
print(f"✅ Saved: {minimized_pdb}")

# Equilibration
print("\n" + "=" * 80)
print("Step 5: Equilibration (100 ps)")
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

simulation.step(50000)  # 100 ps

print("✅ Equilibration complete")

# Production
print("\n" + "=" * 80)
print("Step 6: Production MD (5 ns)")
print("=" * 80)

simulation.reporters.clear()

simulation.reporters.append(
    app.StateDataReporter(
        os.path.join(output_dir, 'production.log'),
        5000,
        step=True,
        time=True,
        potentialEnergy=True,
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

print("  Running 5 ns production...")
simulation.step(2500000)  # 5 ns

# Save final
state = simulation.context.getState(getPositions=True)
final_pdb = os.path.join(output_dir, 'final.pdb')
with open(final_pdb, 'w') as f:
    app.PDBFile.writeFile(modeller.topology, state.getPositions(), f)

print("\n" + "=" * 80)
print("✅ SIMULATION COMPLETE!")
print("=" * 80)
print(f"\n📁 Results: {output_dir}")
print("  - solvated.pdb")
print("  - minimized.pdb")
print("  - equilibration.dcd")
print("  - production.dcd")
print("  - production.log")
print("  - final.pdb")
print("\n" + "=" * 80)
