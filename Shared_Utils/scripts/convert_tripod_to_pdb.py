#!/usr/bin/env python3
"""
Convert Tripod SDF to PDB format
"""

from rdkit import Chem
from rdkit.Chem import AllChem

# Read SDF file
mol = Chem.SDMolSupplier('Tripod.sdf', removeHs=False)[0]

if mol is None:
    print("Error: Could not read SDF file")
    exit(1)

# Add hydrogens if not present
if mol.GetNumAtoms() == mol.GetNumHeavyAtoms():
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)

# Write PDB file
writer = Chem.PDBWriter('tripod_only.pdb')
writer.write(mol)
writer.close()

print(f"Successfully converted Tripod.sdf to tripod_only.pdb")
print(f"Number of atoms: {mol.GetNumAtoms()}")
print(f"Number of heavy atoms: {mol.GetNumHeavyAtoms()}")
print(f"Molecular formula: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
print(f"Molecular weight: {Chem.rdMolDescriptors.CalcExactMolWt(mol):.2f} Da")
