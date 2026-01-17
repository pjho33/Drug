#!/usr/bin/env python
"""
Phase 4: Convert Glycan Residue Names to GLYCAM Format
=======================================================
PDB standard names (NAG, MAN, BMA) -> AMBER GLYCAM names (4YB, 4gA, VMB, 0MA, etc.)

GLYCAM naming convention:
- Position prefix: 0=terminal, 1-6=linkage position
- Sugar code: GA=alpha-GlcNAc, gA=beta-GlcNAc, MA=alpha-Man, mA=beta-Man
- Special: VMB=3,6-branched beta-Man, 4YB=Asn-linked GlcNAc
"""

import sys
import os
from openmm.app import PDBFile

# Paths
INPUT_PDB = "/home/pjho3/projects/Drug/structures/phase4/glut1_glycosylated_man5.pdb"
OUTPUT_PDB = "/home/pjho3/projects/Drug/structures/phase4/glut1_glycam_ready.pdb"


def get_glycam_name(residue, connection_count, is_first_nag):
    """
    Map PDB residue names to GLYCAM names based on connectivity.
    
    GLYCAM naming:
    - 4YB: Asn-linked GlcNAc (first NAG)
    - 4gA: 4-linked beta-GlcNAc
    - VMB: 3,6-branched beta-Mannose (core)
    - 0MA: Terminal alpha-Mannose
    - 3MA: 3-linked alpha-Mannose
    - 6MA: 6-linked alpha-Mannose
    """
    # 1. N-Acetylglucosamine (NAG)
    if residue.name == 'NAG':
        if is_first_nag:
            return '4YB'  # Asn-linked GlcNAc (AMBER specific)
        else:
            return '4gA'  # 4-linked beta-GlcNAc

    # 2. Beta-Mannose (BMA) - usually the core branching point
    if residue.name == 'BMA':
        if connection_count >= 3:  # NAG + Man(α1-3) + Man(α1-6)
            return 'VMB'  # 3,6-branched beta-Mannose
        else:
            return '3mA'  # 3-linked beta-Mannose

    # 3. Alpha-Mannose (MAN)
    if residue.name == 'MAN':
        if connection_count == 1:
            return '0MA'  # Terminal alpha-Mannose
        elif connection_count == 2:
            return '3MA'  # Internal alpha-Mannose (assume 3-linked)
            # Note: Could be 6MA if 6-linked, but requires geometric analysis
            
    return residue.name  # Keep original if no mapping


def find_connections_by_distance(topology, positions, cutoff=0.2):
    """
    Find inter-residue connections based on atomic distances.
    Cutoff in nm (0.2 nm = 2.0 Å, typical covalent bond length)
    """
    import numpy as np
    
    # Get positions as numpy array
    pos_array = np.array([[p.x, p.y, p.z] for p in positions])
    
    # Build residue connection map
    residue_connections = {res: set() for res in topology.residues()}
    
    # Get atoms by residue
    residue_atoms = {}
    for atom in topology.atoms():
        if atom.residue not in residue_atoms:
            residue_atoms[atom.residue] = []
        residue_atoms[atom.residue].append(atom.index)
    
    # Check distances between atoms of different residues
    residue_list = list(topology.residues())
    for i, res1 in enumerate(residue_list):
        for res2 in residue_list[i+1:]:
            # Only check glycan-related connections
            glycan_names = ['NAG', 'MAN', 'BMA', 'ASN']
            if res1.name not in glycan_names and res2.name not in glycan_names:
                continue
            
            # Check if any atoms are within bonding distance
            for idx1 in residue_atoms.get(res1, []):
                for idx2 in residue_atoms.get(res2, []):
                    dist = np.linalg.norm(pos_array[idx1] - pos_array[idx2])
                    if dist < cutoff:
                        residue_connections[res1].add(res2)
                        residue_connections[res2].add(res1)
                        break
                else:
                    continue
                break
    
    return residue_connections


def main():
    print("=" * 60)
    print("Converting Glycan Names to GLYCAM Format")
    print("=" * 60)
    
    # Load PDB
    print(f"\n📂 Loading: {INPUT_PDB}")
    pdb = PDBFile(INPUT_PDB)
    topology = pdb.topology
    positions = pdb.positions
    
    print(f"   Atoms: {topology.getNumAtoms()}")
    
    # Count glycan residues before conversion
    glycan_counts = {}
    for res in topology.residues():
        if res.name in ['NAG', 'MAN', 'BMA']:
            glycan_counts[res.name] = glycan_counts.get(res.name, 0) + 1
    
    print(f"\n📊 Glycan residues found:")
    for name, count in glycan_counts.items():
        print(f"   {name}: {count}")
    
    print("\n🔄 Analyzing connectivity by distance...")
    print("-" * 50)
    
    # Find connections by distance (since CONECT records may be missing)
    residue_connections = find_connections_by_distance(topology, positions, cutoff=0.2)
    
    # Analyze and rename each glycan residue
    nags_found = 0
    conversions = []
    
    for residue in topology.residues():
        # Only process glycan residues
        if residue.name not in ['NAG', 'MAN', 'BMA']:
            continue
        
        # Get connected residues from distance-based analysis
        connected_residues = residue_connections.get(residue, set())
        conn_count = len(connected_residues)
        
        # Check if this is the first NAG (connected to ASN)
        is_first = False
        if residue.name == 'NAG':
            nags_found += 1
            for neighbor_res in connected_residues:
                if neighbor_res.name == 'ASN':
                    is_first = True
                    break
        
        # Get new GLYCAM name
        old_name = residue.name
        new_name = get_glycam_name(residue, conn_count, is_first)
        
        connected_names = [r.name for r in connected_residues]
        print(f"   Residue {residue.index:3d} [{old_name}] -> [{new_name}] "
              f"(Connections: {conn_count}, to: {connected_names})")
        
        if new_name != old_name:
            residue.name = new_name
            conversions.append((old_name, new_name))
    
    print("-" * 50)
    print(f"\n✅ Converted {len(conversions)} residues")
    
    # Save converted PDB
    print(f"\n💾 Saving: {OUTPUT_PDB}")
    with open(OUTPUT_PDB, 'w') as f:
        PDBFile.writeFile(topology, positions, f)
    
    print("\n" + "=" * 60)
    print("Conversion Complete!")
    print("=" * 60)
    print(f"\nOutput file: {OUTPUT_PDB}")
    print("\nNext step: Test with AMBER GLYCAM forcefield")
    print("  ForceField('amber14-all.xml', 'amber14/GLYCAM_06j-1.xml')")


if __name__ == "__main__":
    main()
