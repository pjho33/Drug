#!/usr/bin/env python
import os, sys

# CHARMM-GUI helper 모듈 경로 추가
CHARMM_DIR = '/home/pjho3tr/Downloads/charmm-gui-6750265216membranebuilder/openmm'
sys.path.insert(0, CHARMM_DIR)

from omm_readparams import read_params, read_top
from openmm.app import *
from openmm import *
from openmm.unit import *

SYSTEM_DIR = '/home/pjho3tr/projects/Drug/phase3_glycosylation/glycosylated_new_final'
GPU_INDEX = '0'
CHECKPOINT_FILE = 'prod_glyco_base.chk'
OUTPUT_DCD = 'prod_glyco_base.dcd'
OUTPUT_LOG = 'prod_glyco_base_log.csv'
TOTAL_STEPS = 50000000
TEMPERATURE = 310*kelvin

print("="*80)
print("GLYCOSYLATED BASE SYSTEM - GPU 0")
print("="*80)

os.chdir(SYSTEM_DIR)
resume = os.path.exists(CHECKPOINT_FILE)
print(f"Status: {'RESUMING' if resume else 'NEW SIMULATION'}")

print("\nLoading CHARMM PSF/PDB...")
cwd = os.getcwd()
os.chdir(CHARMM_DIR)
psf = read_top('step5_input.psf', 'CHARMM')
pdb = PDBFile('step5_input.pdb')

print("Loading CHARMM parameters (CHARMM-GUI helper)...")
params = read_params('toppar.str')
os.chdir(cwd)
print("  ✓ Parameters loaded")

print(f"\nPDB atoms: {pdb.topology.getNumAtoms()}")
print("Creating system...")

system = psf.createSystem(
    params,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0*nanometer,
    constraints=HBonds,
    rigidWater=True,
    hydrogenMass=1.5*amu
)
print(f"✓ System created! Particles: {system.getNumParticles()}")

if system.getNumParticles() != pdb.topology.getNumAtoms():
    raise RuntimeError(f"Particle mismatch: {system.getNumParticles()} vs {pdb.topology.getNumAtoms()}")

integrator = LangevinMiddleIntegrator(TEMPERATURE, 1/picosecond, 0.002*picoseconds)
platform = Platform.getPlatformByName('CUDA')
properties = {'DeviceIndex': GPU_INDEX, 'Precision': 'mixed'}
simulation = Simulation(pdb.topology, system, integrator, platform, properties)

def infer_current_step(csv_path):
    if not os.path.exists(csv_path): return 0
    try:
        with open(csv_path, 'r') as f:
            for line in reversed(f.readlines()):
                if line.strip() and not line.startswith('#') and 'Step' not in line:
                    try: return int(line.split(',')[0])
                    except: continue
    except: pass
    return 0

if resume:
    print("\nLoading checkpoint...")
    simulation.loadCheckpoint(CHECKPOINT_FILE)
    current_step = infer_current_step(OUTPUT_LOG)
    remaining_steps = max(0, TOTAL_STEPS - current_step)
    print(f"Resume from: {current_step*0.002/1000:.2f} ns")
else:
    simulation.context.setPositions(pdb.positions)
    print("\nMinimizing...")
    simulation.minimizeEnergy(maxIterations=1000)
    print("Equilibrating 50 ps...")
    simulation.step(25000)
    remaining_steps = TOTAL_STEPS

simulation.reporters.append(DCDReporter(OUTPUT_DCD, 5000, append=resume))
simulation.reporters.append(StateDataReporter(OUTPUT_LOG, 5000,
    step=True, time=True, potentialEnergy=True, temperature=True, speed=True, append=resume))
simulation.reporters.append(CheckpointReporter(CHECKPOINT_FILE, 50000))
simulation.reporters.append(StateDataReporter(sys.stdout, 100000,
    step=True, time=True, speed=True, remainingTime=True, totalSteps=remaining_steps))

print(f"\nRunning {remaining_steps*0.002/1000:.0f} ns production...")
simulation.step(remaining_steps)

final_pdb = 'prod_glyco_base_final.pdb'
with open(final_pdb, 'w') as f:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)
print(f"\n✓ Complete! Saved: {final_pdb}")
