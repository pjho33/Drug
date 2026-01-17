#!/bin/bash
# Phase 3: Glycosylated vs Control with Tripod
# Using Phase 2 validated AMBER force field

source /home/pjho3tr/miniforge3/etc/profile.d/conda.sh
conda activate drug-md

BASE_DIR="/home/pjho3tr/projects/Drug/phase3_glycosylation"

# Phase 2 스크립트 복사
cp /home/pjho3tr/projects/Drug/results/phase2_rep1/gaff_cache.json $BASE_DIR/

# Glycosylated system
echo "=========================================="
echo "Running Glycosylated system..."
echo "=========================================="
cd $BASE_DIR/glycosylated_amber

python << 'PYEOF'
from openmm.app import *
from openmm import *
from openmm.unit import *
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator
import mdtraj as md

# Load protein
pdb = PDBFile('protein.pdb')

# Load ligand
ligand_mol = Molecule.from_file('ligand.sdf')

# Force field with GAFF
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
gaff = GAFFTemplateGenerator(molecules=ligand_mol, cache='../gaff_cache.json')
forcefield.registerTemplateGenerator(gaff.generator)

# Modeller: add ligand and water
modeller = Modeller(pdb.topology, pdb.positions)
modeller.add(ligand_mol.to_topology().to_openmm(), ligand_mol.conformers[0])
modeller.addSolvent(forcefield, model='tip3p', padding=1.0*nanometer)

# Create system
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0*nanometer,
    constraints=HBonds
)

# Simulation
integrator = LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 0.002*picoseconds)
platform = Platform.getPlatformByName('CUDA')
properties = {'DeviceIndex': '0', 'Precision': 'mixed'}

simulation = Simulation(modeller.topology, system, integrator, platform, properties)
simulation.context.setPositions(modeller.positions)

# Minimize
print("Minimizing...")
simulation.minimizeEnergy(maxIterations=1000)

# Equilibrate
print("Equilibrating...")
simulation.step(25000)  # 50ps

# Production
print("Running 100ns production...")
simulation.reporters.append(DCDReporter('prod_glyco_tripod.dcd', 5000))
simulation.reporters.append(StateDataReporter('prod_glyco_tripod_log.csv', 5000,
    step=True, time=True, potentialEnergy=True, temperature=True, speed=True))
simulation.reporters.append(CheckpointReporter('prod_glyco_tripod.chk', 50000))

simulation.step(50000000)  # 100ns

# Save
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('prod_glyco_tripod_final.pdb', 'w'))
print("✓ Glycosylated system complete!")
PYEOF

# Control system
echo "=========================================="
echo "Running Control system..."
echo "=========================================="
cd $BASE_DIR/control_amber

python << 'PYEOF'
from openmm.app import *
from openmm import *
from openmm.unit import *
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator

pdb = PDBFile('protein.pdb')
ligand_mol = Molecule.from_file('ligand.sdf')

forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
gaff = GAFFTemplateGenerator(molecules=ligand_mol, cache='../gaff_cache.json')
forcefield.registerTemplateGenerator(gaff.generator)

modeller = Modeller(pdb.topology, pdb.positions)
modeller.add(ligand_mol.to_topology().to_openmm(), ligand_mol.conformers[0])
modeller.addSolvent(forcefield, model='tip3p', padding=1.0*nanometer)

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
print("✓ Control system complete!")
PYEOF

echo "=========================================="
echo "Phase 3 simulations complete!"
echo "=========================================="
