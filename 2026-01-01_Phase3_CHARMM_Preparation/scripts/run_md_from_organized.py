#!/usr/bin/env python
import sys
sys.path.insert(0, '/home/pjho3tr/projects/Drug/phase3_glycosylation/glycosylated')

from openmm.app import *
from openmm import *
from openmm.unit import *
import omm_readinputs
import omm_readparams

# Load system
print("Loading system...")
psf = CharmmPsfFile('system.psf')
pdb = PDBFile('complex_with_tripod.pdb')

# Load parameters
print("Loading CHARMM parameters...")
params = CharmmParameterSet('toppar.str')

# Create system
print("Creating system...")
system = psf.createSystem(
    params,
    nonbondedMethod=PME,
    nonbondedCutoff=1.2*nanometer,
    switchDistance=1.0*nanometer,
    constraints=HBonds
)

# Integrator
integrator = LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 0.002*picoseconds)

# Platform
try:
    platform = Platform.getPlatformByName('CUDA')
    properties = {'DeviceIndex': '0', 'Precision': 'mixed'}
    print("Using CUDA")
except:
    platform = Platform.getPlatformByName('CPU')
    properties = {}
    print("Using CPU")

# Simulation
simulation = Simulation(psf.topology, system, integrator, platform, properties)
simulation.context.setPositions(pdb.positions)

# Minimize
print("Minimizing...")
simulation.minimizeEnergy(maxIterations=1000)

# Reporters
simulation.reporters.append(DCDReporter('prod_glycosylated_tripod.dcd', 5000))
simulation.reporters.append(StateDataReporter(
    'prod_glycosylated_tripod_log.csv', 5000,
    step=True, time=True, potentialEnergy=True, temperature=True, speed=True
))
simulation.reporters.append(CheckpointReporter('prod_glycosylated_tripod.chk', 50000))

# Run 100ns
print("Running 100ns MD...")
simulation.step(50000000)

# Save
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('prod_glycosylated_tripod_final.pdb', 'w'))
print("Complete!")
