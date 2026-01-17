#!/usr/bin/env python3
"""
Control Group GLUT1-SDG Complex MD Simulation
Using OpenMM with CHARMM force field
"""

import openmm as mm
from openmm import app, unit
import sys
import os
from datetime import datetime

print("=" * 80)
print("Control Group: GLUT1-SDG Complex MD Simulation")
print("=" * 80)
print(f"Start time: {datetime.now()}")

# Check CUDA
try:
    platform = mm.Platform.getPlatformByName('CUDA')
    print(f"✅ CUDA available")
    use_cuda = True
except:
    platform = mm.Platform.getPlatformByName('CPU')
    print(f"⚠️  Using CPU (slower)")
    use_cuda = False

# File paths
work_dir = '/home/pjho3/projects/Drug/control_md_simulation'
psf_file = f'{work_dir}/step4_lipid.psf'
pdb_file = f'{work_dir}/step4_lipid.pdb'
toppar_str = f'{work_dir}/toppar.str'
output_dir = f'{work_dir}/results'

os.makedirs(output_dir, exist_ok=True)

print(f"\n📁 Working directory: {work_dir}")
print(f"📁 Output directory: {output_dir}")

# Step 1: Load PSF and PDB
print("\n" + "=" * 80)
print("Step 1: Loading system")
print("=" * 80)

print(f"  Loading PSF: {psf_file}")
psf = app.CharmmPsfFile(psf_file)
print(f"  ✅ PSF loaded: {psf.topology.getNumAtoms()} atoms")

print(f"  Loading PDB: {pdb_file}")
pdb = app.PDBFile(pdb_file)
print(f"  ✅ PDB loaded")

# Count components
protein_atoms = 0
lipid_atoms = 0
water_atoms = 0
ion_atoms = 0
sdg_atoms = 0

for residue in psf.topology.residues():
    res_name = residue.name
    n_atoms = len(list(residue.atoms()))
    
    if res_name == 'SDG':
        sdg_atoms += n_atoms
    elif res_name in ['TIP3', 'HOH']:
        water_atoms += n_atoms
    elif res_name in ['POT', 'CLA', 'SOD', 'CL']:
        ion_atoms += n_atoms
    elif res_name in ['POPC', 'POPE', 'POPS', 'CHOL']:
        lipid_atoms += n_atoms
    else:
        protein_atoms += n_atoms

print(f"\n  System composition:")
print(f"    Protein: {protein_atoms} atoms")
print(f"    SDG ligand: {sdg_atoms} atoms")
print(f"    Lipids: {lipid_atoms} atoms")
print(f"    Water: {water_atoms} atoms")
print(f"    Ions: {ion_atoms} atoms")

# Step 2: Load parameters
print("\n" + "=" * 80)
print("Step 2: Loading CHARMM parameters")
print("=" * 80)

toppar_dir = f'{work_dir}/toppar'
print(f"  Loading individual parameter files from: {toppar_dir}")

# Load parameters in correct order
params = app.CharmmParameterSet(
    f'{toppar_dir}/top_all36_prot.rtf',
    f'{toppar_dir}/top_all36_lipid.rtf',
    f'{toppar_dir}/top_all36_carb.rtf',
    f'{toppar_dir}/top_all36_cgenff.rtf',
    f'{toppar_dir}/sdg.rtf',
    f'{toppar_dir}/par_all36m_prot.prm',
    f'{toppar_dir}/par_all36_lipid.prm',
    f'{toppar_dir}/par_all36_carb.prm',
    f'{toppar_dir}/par_all36_cgenff.prm',
    f'{toppar_dir}/sdg.prm'
)
print(f"  ✅ Parameters loaded (including SDG)")

# Step 3: Create system
print("\n" + "=" * 80)
print("Step 3: Creating system")
print("=" * 80)

system = psf.createSystem(
    params,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.2*unit.nanometer,
    constraints=app.HBonds,
    rigidWater=True,
    ewaldErrorTolerance=0.0005
)

print(f"  ✅ System created")
print(f"  Forces: {system.getNumForces()}")

# Step 4: Energy minimization
print("\n" + "=" * 80)
print("Step 4: Energy minimization")
print("=" * 80)

integrator = mm.LangevinMiddleIntegrator(
    310*unit.kelvin,
    1.0/unit.picosecond,
    2.0*unit.femtosecond
)

if use_cuda:
    properties = {'CudaPrecision': 'mixed', 'DeviceIndex': '0'}
else:
    properties = {}

simulation = app.Simulation(
    psf.topology,
    system,
    integrator,
    platform,
    properties
)

simulation.context.setPositions(pdb.positions)

print("  Initial energy...")
state = simulation.context.getState(getEnergy=True)
initial_energy = state.getPotentialEnergy()
print(f"    Potential Energy: {initial_energy}")

print("  Minimizing (max 5000 steps)...")
simulation.minimizeEnergy(maxIterations=5000)

state = simulation.context.getState(getEnergy=True, getPositions=True)
final_energy = state.getPotentialEnergy()
print(f"    Final Potential Energy: {final_energy}")
print(f"    Energy change: {final_energy - initial_energy}")

# Save minimized structure
minimized_pdb = f'{output_dir}/minimized.pdb'
with open(minimized_pdb, 'w') as f:
    app.PDBFile.writeFile(psf.topology, state.getPositions(), f)
print(f"  ✅ Saved: {minimized_pdb}")

# Step 5: Equilibration (NVT)
print("\n" + "=" * 80)
print("Step 5: Equilibration (NVT, 500 ps)")
print("=" * 80)

simulation.context.setVelocitiesToTemperature(310*unit.kelvin)

simulation.reporters.append(
    app.StateDataReporter(
        sys.stdout,
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
    app.DCDReporter(f'{output_dir}/equilibration.dcd', 5000)
)

simulation.reporters.append(
    app.StateDataReporter(
        f'{output_dir}/equilibration.log',
        1000,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        temperature=True,
        volume=True,
        density=True
    )
)

print("  Running 500 ps equilibration...")
simulation.step(250000)  # 500 ps at 2 fs timestep

print("  ✅ Equilibration complete")

# Save equilibrated structure
equilibrated_pdb = f'{output_dir}/equilibrated.pdb'
state = simulation.context.getState(getPositions=True)
with open(equilibrated_pdb, 'w') as f:
    app.PDBFile.writeFile(psf.topology, state.getPositions(), f)
print(f"  ✅ Saved: {equilibrated_pdb}")

# Step 6: Production MD
print("\n" + "=" * 80)
print("Step 6: Production MD (10 ns)")
print("=" * 80)

simulation.reporters.clear()

simulation.reporters.append(
    app.StateDataReporter(
        sys.stdout,
        10000,
        step=True,
        time=True,
        potentialEnergy=True,
        temperature=True,
        speed=True,
        remainingTime=True,
        totalSteps=5000000
    )
)

simulation.reporters.append(
    app.DCDReporter(f'{output_dir}/production.dcd', 5000)  # Save every 10 ps
)

simulation.reporters.append(
    app.StateDataReporter(
        f'{output_dir}/production.log',
        5000,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=True,
        density=True
    )
)

simulation.reporters.append(
    app.CheckpointReporter(f'{output_dir}/checkpoint.chk', 50000)
)

print("  Running 10 ns production...")
print("  (This will take several hours...)")
simulation.step(5000000)  # 10 ns at 2 fs timestep

# Save final structure
final_pdb = f'{output_dir}/final.pdb'
state = simulation.context.getState(getPositions=True)
with open(final_pdb, 'w') as f:
    app.PDBFile.writeFile(psf.topology, state.getPositions(), f)

print("\n" + "=" * 80)
print("✅ SIMULATION COMPLETE!")
print("=" * 80)
print(f"End time: {datetime.now()}")
print(f"\n📁 Results saved in: {output_dir}")
print("\nOutput files:")
print(f"  - minimized.pdb")
print(f"  - equilibration.dcd")
print(f"  - equilibration.log")
print(f"  - equilibrated.pdb")
print(f"  - production.dcd")
print(f"  - production.log")
print(f"  - checkpoint.chk")
print(f"  - final.pdb")
print("\n" + "=" * 80)
