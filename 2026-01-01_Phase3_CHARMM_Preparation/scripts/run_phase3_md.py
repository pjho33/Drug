#!/usr/bin/env python
"""
Phase 3: Glycosylation Effect MD Simulation
Tripod binding to Glycosylated vs Control GLUT1
"""
from openmm.app import *
from openmm import *
from openmm.unit import *
import sys
import os

def run_md_simulation(system_name, pdb_file, output_prefix, steps=50000000):
    """
    Run MD simulation for Phase 3
    system_name: 'glycosylated' or 'control'
    steps: 50M steps = 100ns (2fs timestep)
    """
    print(f"\n{'='*80}")
    print(f"Starting MD simulation: {system_name.upper()}")
    print(f"{'='*80}")
    
    # Load PDB
    print(f"Loading {pdb_file}...")
    pdb = PDBFile(pdb_file)
    
    # Force field
    print("Setting up force field...")
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    
    # Create system
    print("Creating system...")
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0*nanometer,
        constraints=HBonds,
        rigidWater=True,
        ewaldErrorTolerance=0.0005
    )
    
    # Integrator
    integrator = LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 0.002*picoseconds)
    
    # Platform
    try:
        platform = Platform.getPlatformByName('CUDA')
        properties = {'DeviceIndex': '0', 'Precision': 'mixed'}
        print("Using CUDA platform")
    except:
        try:
            platform = Platform.getPlatformByName('OpenCL')
            properties = {'DeviceIndex': '0', 'Precision': 'mixed'}
            print("Using OpenCL platform")
        except:
            platform = Platform.getPlatformByName('CPU')
            properties = {}
            print("Using CPU platform")
    
    # Simulation
    simulation = Simulation(pdb.topology, system, integrator, platform, properties)
    simulation.context.setPositions(pdb.positions)
    
    # Minimize
    print("Minimizing energy...")
    simulation.minimizeEnergy(maxIterations=1000)
    
    # Reporters
    simulation.reporters.append(DCDReporter(f'{output_prefix}.dcd', 5000))
    simulation.reporters.append(StateDataReporter(
        f'{output_prefix}_log.csv', 5000,
        step=True, time=True, potentialEnergy=True, kineticEnergy=True,
        totalEnergy=True, temperature=True, volume=True, density=True,
        speed=True, separator=','
    ))
    simulation.reporters.append(CheckpointReporter(f'{output_prefix}.chk', 50000))
    
    # Run
    print(f"Running {steps} steps ({steps*0.002/1000000:.1f} ns)...")
    simulation.step(steps)
    
    # Save final state
    print("Saving final structure...")
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, open(f'{output_prefix}_final.pdb', 'w'))
    
    print(f"✓ Simulation complete: {system_name}")

if __name__ == "__main__":
    base_dir = "/home/pjho3tr/projects/Drug/phase3_glycosylation"
    
    # Glycosylated system
    print("\n" + "="*80)
    print("PHASE 3: GLYCOSYLATION EFFECT ON TRIPOD BINDING")
    print("="*80)
    
    systems = [
        ('glycosylated', f'{base_dir}/glycosylated/complex_with_tripod.pdb', 
         f'{base_dir}/glycosylated/prod_glyco_tripod'),
        ('control', f'{base_dir}/control/complex_with_tripod.pdb',
         f'{base_dir}/control/prod_control_tripod')
    ]
    
    for system_name, pdb_file, output_prefix in systems:
        if os.path.exists(pdb_file):
            run_md_simulation(system_name, pdb_file, output_prefix, steps=50000000)
        else:
            print(f"ERROR: {pdb_file} not found!")
    
    print("\n" + "="*80)
    print("All simulations complete!")
    print("="*80)
