#!/usr/bin/env python
"""
Phase 4: Load CHARMM-GUI Generated System for OpenMM
=====================================================
Load the glycosylated GLUT1 system from CHARMM-GUI output
and prepare it for OpenMM simulations.
"""

import os
from openmm.app import *
from openmm import *
from openmm.unit import *
import mdtraj as md

# Paths
CHARMM_GUI_DIR = "/home/pjho3/projects/Drug/structures/phase4/charmm_gui_output"
GROMACS_DIR = os.path.join(CHARMM_GUI_DIR, "gromacs")
OUTPUT_DIR = "/home/pjho3/projects/Drug/structures/phase4"

def load_charmm_gui_system():
    """Load CHARMM-GUI generated system"""
    
    print("=" * 60)
    print("Loading CHARMM-GUI System for OpenMM")
    print("=" * 60)
    
    # Load structure files
    gro_file = os.path.join(GROMACS_DIR, "step5_input.gro")
    top_file = os.path.join(GROMACS_DIR, "topol.top")
    
    print(f"\n📂 Loading structure...")
    print(f"   GRO: {gro_file}")
    print(f"   TOP: {top_file}")
    
    # Check files exist
    if not os.path.exists(gro_file):
        raise FileNotFoundError(f"GRO file not found: {gro_file}")
    if not os.path.exists(top_file):
        raise FileNotFoundError(f"TOP file not found: {top_file}")
    
    # Load with OpenMM
    print("\n🔧 Loading with OpenMM GromacsGroFile...")
    gro = GromacsGroFile(gro_file)
    
    print(f"   Atoms: {len(gro.positions)}")
    print(f"   Box: {gro.getPeriodicBoxVectors()}")
    
    # Load topology
    print("\n🔧 Loading topology...")
    top = GromacsTopFile(top_file, 
                         periodicBoxVectors=gro.getPeriodicBoxVectors(),
                         includeDir=os.path.join(GROMACS_DIR, "toppar"))
    
    print(f"   Residues: {top.topology.getNumResidues()}")
    print(f"   Atoms: {top.topology.getNumAtoms()}")
    print(f"   Bonds: {top.topology.getNumBonds()}")
    
    # Analyze composition
    print("\n📊 System composition:")
    residue_counts = {}
    for res in top.topology.residues():
        if res.name not in residue_counts:
            residue_counts[res.name] = 0
        residue_counts[res.name] += 1
    
    # Show key components
    key_residues = ['PROA', 'TIP3', 'POPC', 'POPE', 'POPS', 'CHL1', 'PSM', 'POT', 'CLA']
    for name in key_residues:
        if name in residue_counts:
            print(f"   {name}: {residue_counts[name]}")
    
    # Check for glycan residues
    glycan_names = [name for name in residue_counts.keys() 
                    if any(x in name for x in ['NAG', 'MAN', 'BMA', 'AMAN', 'BMAN', 'BGNA'])]
    if glycan_names:
        print(f"\n   Glycan residues found:")
        for name in glycan_names:
            print(f"      {name}: {residue_counts[name]}")
    
    # Save PDB for visualization
    print("\n💾 Saving PDB for visualization...")
    pdb_out = os.path.join(OUTPUT_DIR, "glut1_glycosylated_charmm_gui.pdb")
    PDBFile.writeFile(top.topology, gro.positions, open(pdb_out, 'w'))
    print(f"   Saved: {pdb_out}")
    
    return top, gro


def create_openmm_system(top, gro):
    """Create OpenMM system from CHARMM-GUI topology"""
    
    print("\n" + "=" * 60)
    print("Creating OpenMM System")
    print("=" * 60)
    
    print("\n🔧 Creating system...")
    system = top.createSystem(nonbondedMethod=PME,
                              nonbondedCutoff=1.2*nanometer,
                              constraints=HBonds,
                              rigidWater=True,
                              ewaldErrorTolerance=0.0005)
    
    print(f"✅ System created!")
    print(f"   Particles: {system.getNumParticles()}")
    print(f"   Forces: {system.getNumForces()}")
    
    # Show forces
    print("\n   Force types:")
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        print(f"      {i}: {force.__class__.__name__}")
    
    return system


def test_simulation(system, top, gro):
    """Test a short simulation"""
    
    print("\n" + "=" * 60)
    print("Testing Short Simulation")
    print("=" * 60)
    
    # Setup integrator
    temperature = 310*kelvin
    friction = 1/picosecond
    timestep = 2*femtoseconds
    
    integrator = LangevinMiddleIntegrator(temperature, friction, timestep)
    
    # Create simulation
    print("\n🔧 Creating simulation...")
    platform = Platform.getPlatformByName('CPU')
    simulation = Simulation(top.topology, system, integrator, platform)
    simulation.context.setPositions(gro.positions)
    
    # Energy minimization
    print("\n⚡ Energy minimization...")
    print(f"   Initial energy: {simulation.context.getState(getEnergy=True).getPotentialEnergy()}")
    
    simulation.minimizeEnergy(maxIterations=100)
    
    print(f"   Final energy: {simulation.context.getState(getEnergy=True).getPotentialEnergy()}")
    
    # Short test run
    print("\n🏃 Test run (100 steps)...")
    simulation.step(100)
    
    print("✅ Test simulation successful!")
    
    return simulation


def main():
    """Main function"""
    
    try:
        # Load system
        top, gro = load_charmm_gui_system()
        
        # Create OpenMM system
        system = create_openmm_system(top, gro)
        
        # Test simulation
        simulation = test_simulation(system, top, gro)
        
        print("\n" + "=" * 60)
        print("✨ CHARMM-GUI System Successfully Loaded!")
        print("=" * 60)
        print("\nNext step: Set up penetration test simulations")
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
