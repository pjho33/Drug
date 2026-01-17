#!/usr/bin/env python
"""
Phase 3: Glycosylated GLUT1 + Tripod MD Simulation (GPU 0)
PSF + CRD 방식 사용 (PDB 포맷 제한 회피)
"""
import os, sys
sys.path.insert(0, '/home/pjho3tr/Downloads/charmm-gui-6750265216membranebuilder/openmm')

from omm_readparams import read_params, read_top, read_crd
from openmm.app import *
from openmm import *
from openmm.unit import *

CHARMM_DIR = '/home/pjho3tr/Downloads/charmm-gui-6750265216membranebuilder/openmm'
WORK_DIR = '/home/pjho3tr/projects/Drug/phase3_with_tripod/experimental'
GPU_INDEX = '0'

TOTAL_STEPS = 50000000
TEMPERATURE = 310*kelvin
FRICTION = 1/picosecond
DT = 0.002*picoseconds

CHECKPOINT_FILE = 'prod_experimental.chk'
OUTPUT_DCD = 'prod_experimental.dcd'
OUTPUT_LOG = 'prod_experimental.log'
FINAL_PDB = 'prod_experimental_final.pdb'

print("="*80)
print("EXPERIMENTAL: Glycosylated GLUT1 + Tripod (GPU 0)")
print("="*80)

os.chdir(WORK_DIR)
resume = os.path.exists(CHECKPOINT_FILE)
print(f"Status: {'RESUMING' if resume else 'NEW SIMULATION'}")

print("\nLoading system...")
os.chdir(CHARMM_DIR)

# PSF + CRD 사용 (PDB 대신)
psf = read_top('step5_input.psf', 'CHARMM')
crd = read_crd('step5_input.crd', 'CHARMM')

print("Loading parameters (including Tripod)...")
params = read_params('toppar.str')
os.chdir(WORK_DIR)

print(f"  Topology atoms: {psf.topology.getNumAtoms()}")

# Tripod 좌표를 CRD에 추가해야 함
# 간단한 방법: 원본 PDB 사용 (Tripod 포함)
print("\nLoading Tripod coordinates from PDB...")
tripod_pdb = PDBFile('/home/pjho3tr/projects/Drug/phase3_with_tripod/experimental/step5_input_with_tripod_fixed.pdb')

print("\nCreating system...")
system = psf.createSystem(
    params,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0*nanometer,
    constraints=HBonds,
    rigidWater=True,
    hydrogenMass=1.5*amu
)
print(f"  ✓ System particles: {system.getNumParticles()}")

integrator = LangevinMiddleIntegrator(TEMPERATURE, FRICTION, DT)
platform = Platform.getPlatformByName('CUDA')
properties = {'DeviceIndex': GPU_INDEX, 'Precision': 'mixed'}

simulation = Simulation(tripod_pdb.topology, system, integrator, platform, properties)

if resume:
    print("\nLoading checkpoint...")
    simulation.loadCheckpoint(CHECKPOINT_FILE)
    current_step = 0
    if os.path.exists(OUTPUT_LOG):
        with open(OUTPUT_LOG, 'r') as f:
            for line in reversed(f.readlines()):
                if line.strip() and not line.startswith('#'):
                    try:
                        current_step = int(line.split(',')[0])
                        break
                    except: continue
    remaining_steps = max(0, TOTAL_STEPS - current_step)
    print(f"  Resume from: {current_step*0.002/1000:.2f} ns")
else:
    print("\nInitializing...")
    simulation.context.setPositions(tripod_pdb.positions)
    print("  Minimizing...")
    simulation.minimizeEnergy(maxIterations=1000)
    print("  Equilibrating (50 ps)...")
    simulation.step(25000)
    remaining_steps = TOTAL_STEPS

simulation.reporters.append(DCDReporter(OUTPUT_DCD, 5000, append=resume))
simulation.reporters.append(StateDataReporter(OUTPUT_LOG, 5000,
    step=True, time=True, potentialEnergy=True, temperature=True, speed=True, append=resume))
simulation.reporters.append(CheckpointReporter(CHECKPOINT_FILE, 50000))
simulation.reporters.append(StateDataReporter(sys.stdout, 100000,
    step=True, time=True, speed=True, remainingTime=True, totalSteps=remaining_steps))

print(f"\nRunning {remaining_steps*0.002/1000:.1f} ns production...")
simulation.step(remaining_steps)

print("\nSaving final structure...")
with open(FINAL_PDB, 'w') as f:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)

print(f"\n✓ COMPLETE! Final: {FINAL_PDB}")
