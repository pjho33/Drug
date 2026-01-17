#!/usr/bin/env python3
"""
Create PSF file for glut1_tripod_complex.pdb using Python
This is a simplified PSF generator for OpenMM compatibility
"""

from openmm.app import PDBFile
import sys

print("=" * 80)
print("Creating PSF for GLUT1-SDG Complex")
print("=" * 80)

# Load PDB
pdb_file = 'glut1_tripod_complex.pdb'
print(f"\n1. Loading {pdb_file}...")

pdb = PDBFile(pdb_file)
topology = pdb.topology

print(f"   ✅ Loaded: {topology.getNumAtoms()} atoms")
print(f"   Residues: {topology.getNumResidues()}")
print(f"   Chains: {topology.getNumChains()}")

# Count residue types
residue_counts = {}
for residue in topology.residues():
    name = residue.name
    residue_counts[name] = residue_counts.get(name, 0) + 1

print(f"\n2. Residue composition:")
for name, count in sorted(residue_counts.items()):
    print(f"   {name}: {count}")

# Write PSF file
psf_file = 'glut1_tripod_complex.psf'
print(f"\n3. Writing PSF file: {psf_file}")

with open(psf_file, 'w') as f:
    # Header
    f.write("PSF CMAP CHEQ\n\n")
    f.write("       1 !NTITLE\n")
    f.write(" REMARKS Generated for GLUT1-SDG complex\n\n")
    
    # Atoms section
    atoms = list(topology.atoms())
    f.write(f"{len(atoms):>8} !NATOM\n")
    
    for i, atom in enumerate(atoms, 1):
        residue = atom.residue
        chain = residue.chain
        
        # Segment ID (chain ID)
        segid = chain.id if chain.id else 'PROA'
        
        # Residue info
        resid = residue.index + 1
        resname = residue.name
        
        # Atom info
        atomname = atom.name
        atomtype = atom.element.symbol if atom.element else 'X'
        
        # Charge and mass (defaults)
        charge = 0.0
        mass = atom.element.mass._value if atom.element else 0.0
        
        # Write atom line
        # Format: atom_id segid resid resname atomname atomtype charge mass
        f.write(f"{i:>8} {segid:<4} {resid:<4} {resname:<4} {atomname:<4} {atomtype:<4} {charge:>10.6f} {mass:>13.4f}           0\n")
    
    # Bonds section
    bonds = list(topology.bonds())
    f.write(f"\n{len(bonds):>8} !NBOND: bonds\n")
    
    bond_count = 0
    for bond in bonds:
        atom1_idx = bond[0].index + 1
        atom2_idx = bond[1].index + 1
        f.write(f"{atom1_idx:>8}{atom2_idx:>8}")
        bond_count += 1
        if bond_count % 4 == 0:
            f.write("\n")
    
    if bond_count % 4 != 0:
        f.write("\n")
    
    # Empty sections (angles, dihedrals, impropers, donors, acceptors, NNB)
    for section in ["!NTHETA: angles", "!NPHI: dihedrals", "!NIMPHI: impropers",
                    "!NDON: donors", "!NACC: acceptors", "!NNB"]:
        f.write(f"\n       0 {section}\n")
    
    # Groups section
    f.write(f"\n       1       0 !NGRP\n")
    f.write(f"       0       0       0\n")

print(f"   ✅ PSF written: {len(atoms)} atoms, {len(bonds)} bonds")

# Verify PSF
print(f"\n4. Verifying PSF...")
try:
    with open(psf_file, 'r') as f:
        lines = f.readlines()
    
    # Check for SDG
    sdg_lines = [l for l in lines if ' SDG ' in l]
    print(f"   ✅ SDG atoms in PSF: {len(sdg_lines)}")
    
    if len(sdg_lines) == 0:
        print("   ⚠️  WARNING: No SDG atoms found!")
    
except Exception as e:
    print(f"   ❌ Error: {e}")

print("\n" + "=" * 80)
print("✅ PSF Creation Complete!")
print("=" * 80)
print(f"\nOutput: {psf_file}")
print("\nNext: Use CharmmPsfFile to load this PSF in OpenMM")
print("=" * 80)
