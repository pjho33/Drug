#!/bin/bash
# Phase 3: Hybrid Force Field (CHARMM36 + GAFF)
source /home/pjho3tr/miniforge3/etc/profile.d/conda.sh
conda activate drug-md

BASE_DIR="/home/pjho3tr/projects/Drug/phase3_glycosylation"
cp /home/pjho3tr/projects/Drug/results/phase2_rep1/gaff_cache.json $BASE_DIR/

echo "=========================================="
echo "Glycosylated System (CHARMM36 + GAFF)"
echo "=========================================="
cd $BASE_DIR/glycosylated_final

python << 'PYEOF'
from openmm.app import *
from openmm import *
from openmm.unit import *
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator

print("Loading protein+ligand complex...")
pdb = PDBFile('system_with_tripod.pdb')

print("Loading ligand for GAFF...")
ligand_mol = Molecule.from_file('ligand.sdf')

print("Setting up CHARMM36 + GAFF hybrid force field...")
forcefield = ForceField('charmm36.xml', 'charmm36/water.xml')
gaff = GAFFTemplateGenerator(molecules=ligand_mol, cache='../gaff_cache.json')
forcefield.registerTemplateGenerator(gaff.generator)

print("Adding solvent...")
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(forcefield, model='tip3p', padding=1.0*nanometer, ionicStrength=0.15*molar)

print("Creating system...")
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0*nanometer,
    constraints=HBonds
)

integrator = LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 0.002*picoseconds)
platform = Platform.getPlatformByName('CUDA')
properties = {'DeviceIndex': '0', 'Precision': 'mixed'}

simulation = Simulation(modeller.topology, system, integrator, platform, properties)
simulation.context.setPositions(modeller.positions)

print("Minimizing...")
simulation.minimizeEnergy(maxIterations=1000)

print("Equilibrating 50ps...")
simulation.step(25000)

print("Running 100ns production...")
simulation.reporters.append(DCDReporter('prod_glyco_tripod.dcd', 5000))
simulation.reporters.append(StateDataReporter('prod_glyco_tripod_log.csv', 5000,
    step=True, time=True, potentialEnergy=True, temperature=True, speed=True))
simulation.reporters.append(CheckpointReporter('prod_glyco_tripod.chk', 50000))

simulation.step(50000000)

positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('prod_glyco_tripod_final.pdb', 'w'))
print("✓ Glycosylated complete!")
PYEOF

echo "=========================================="
echo "Control System (CHARMM36 + GAFF)"
echo "=========================================="
cd $BASE_DIR/control_final

python << 'PYEOF'
from openmm.app import *
from openmm import *
from openmm.unit import *
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator

pdb = PDBFile('system_with_tripod.pdb')
ligand_mol = Molecule.from_file('ligand.sdf')

forcefield = ForceField('charmm36.xml', 'charmm36/water.xml')
gaff = GAFFTemplateGenerator(molecules=ligand_mol, cache='../gaff_cache.json')
forcefield.registerTemplateGenerator(gaff.generator)

modeller = Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(forcefield, model='tip3p', padding=1.0*nanometer, ionicStrength=0.15*molar)

system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME,
    nonbondedCutoff=1.0*nanometer, constraints=HBonds)

integrator = LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 0.002*picoseconds)
platform = Platform.getPlatformByName('CUDA')
properties = {'DeviceIndex': '0', 'Precision': 'mixed'}

simulation = Simulation(modeller.topology, system, integrator, platform, properties)
simulation.context.setPositions(modeller.positions)

print("Minimizing...")
simulation.minimizeEnergy(maxIterations=1000)

print("Equilibrating...")
simulation.step(25000)

print("Running 100ns production...")
simulation.reporters.append(DCDReporter('prod_control_tripod.dcd', 5000))
simulation.reporters.append(StateDataReporter('prod_control_tripod_log.csv', 5000,
    step=True, time=True, potentialEnergy=True, temperature=True, speed=True))
simulation.reporters.append(CheckpointReporter('prod_control_tripod.chk', 50000))

simulation.step(50000000)

positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('prod_control_tripod_final.pdb', 'w'))
print("✓ Control complete!")
PYEOF

echo "=========================================="
echo "Phase 3 Complete!"
echo "=========================================="
