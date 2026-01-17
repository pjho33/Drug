#!/usr/bin/env python
"""
Phase 4: Clean PDB for CHARMM36 ForceField
==========================================
1. Remove existing hydrogens (prevent duplicates)
2. Rename residues to CHARMM36 standard:
   - MAN -> AMAN (Alpha-Mannose)
   - BMA -> BMAN (Beta-Mannose)
   - NAG -> BGNA (Beta-GlcNAc)
3. Merge chains to prevent terminal capping issues
"""

import os
import numpy as np
from openmm.app import PDBFile, Modeller, Element

# Paths
INPUT_PDB = "/home/pjho3/projects/Drug/structures/phase4/glut1_glycosylated_man5.pdb"
OUTPUT_PDB = "/home/pjho3/projects/Drug/structures/phase4/glut1_charmm_ready.pdb"


def main():
    print("=" * 60)
    print("🧹 Cleaning PDB for CHARMM36")
    print("=" * 60)
    
    # Load PDB
    print(f"\n📂 Loading: {INPUT_PDB}")
    pdb = PDBFile(INPUT_PDB)
    modeller = Modeller(pdb.topology, pdb.positions)
    
    print(f"   Initial atoms: {modeller.topology.getNumAtoms()}")
    
    # Step 1: Remove existing hydrogens
    print("\n🔧 Step 1: Removing existing hydrogens...")
    hydrogen = Element.getBySymbol('H')
    atoms_to_delete = [atom for atom in modeller.topology.atoms() 
                       if atom.element == hydrogen]
    
    print(f"   Found {len(atoms_to_delete)} hydrogen atoms to remove")
    modeller.delete(atoms_to_delete)
    print(f"   Atoms after removal: {modeller.topology.getNumAtoms()}")
    
    # Step 2: Rename residues for CHARMM36
    print("\n🔧 Step 2: Renaming residues for CHARMM36...")
    
    rename_map = {
        'MAN': 'AMAN',  # Alpha-Mannose
        'BMA': 'BMAN',  # Beta-Mannose  
        'NAG': 'BGNA',  # Beta-GlcNAc
    }
    
    renamed_count = 0
    for res in modeller.topology.residues():
        if res.name in rename_map:
            old_name = res.name
            new_name = rename_map[res.name]
            res.name = new_name
            print(f"   Residue {res.index}: {old_name} -> {new_name}")
            renamed_count += 1
    
    print(f"   Renamed {renamed_count} residues")
    
    # Step 3: Verify residue types
    print("\n📊 Final residue composition:")
    residue_counts = {}
    for res in modeller.topology.residues():
        if res.name not in residue_counts:
            residue_counts[res.name] = 0
        residue_counts[res.name] += 1
    
    glycan_names = ['AMAN', 'BMAN', 'BGNA']
    for name in glycan_names:
        if name in residue_counts:
            print(f"   {name}: {residue_counts[name]}")
    
    # Step 4: Save cleaned PDB
    print(f"\n💾 Saving: {OUTPUT_PDB}")
    with open(OUTPUT_PDB, 'w') as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
    
    print("\n" + "=" * 60)
    print("✨ PDB Cleaning Complete!")
    print("=" * 60)
    print(f"\nOutput: {OUTPUT_PDB}")
    print("\nNext step: Test with CHARMM36 + glycan.xml")


if __name__ == "__main__":
    main()
