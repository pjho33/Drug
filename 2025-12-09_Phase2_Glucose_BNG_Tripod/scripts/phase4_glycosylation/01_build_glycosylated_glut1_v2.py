#!/usr/bin/env python
"""
Phase 4: Build Glycosylated GLUT1 Model (Template-based)
=========================================================
Uses PDB template from 1GYA (Man5 N-glycan) to attach to GLUT1 Asn45

Workflow:
1. Load Man5 glycan template from 1GYA
2. Load GLUT1 and find Asn45 (THR45 in crystal → restore to ASN)
3. Superposition: Align glycan to Asn45 ND2 position
4. Rotate glycan to point outward from protein surface
5. Merge structures
6. Create bond between ASN45-ND2 and NAG-C1

Glycan structure (from 1GYA):
    NAG-1 (GlcNAc) → NAG-2 (GlcNAc) → BMA-3 (β-Man) → MAN branches
    Total: ~200 atoms, 7 sugar residues
"""

import os
import numpy as np
from scipy.spatial.transform import Rotation
import warnings
warnings.filterwarnings('ignore')

# Paths
RAW_PDB = "/home/pjho3/projects/Drug/raw_data/4PYP.pdb"
GLYCAN_TEMPLATE = "/home/pjho3/projects/Drug/structures/phase4/man5_glycan_template.pdb"
OUTPUT_DIR = "/home/pjho3/projects/Drug/results/phase4_glycan"
STRUCTURES_DIR = "/home/pjho3/projects/Drug/structures/phase4"

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(STRUCTURES_DIR, exist_ok=True)


def parse_pdb(pdb_path):
    """Parse PDB file and return atoms as list of dicts"""
    atoms = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom = {
                    'record': line[:6].strip(),
                    'serial': int(line[6:11]),
                    'name': line[12:16].strip(),
                    'altloc': line[16],
                    'resname': line[17:20].strip(),
                    'chain': line[21],
                    'resseq': int(line[22:26]),
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                    'occupancy': float(line[54:60]) if len(line) > 54 else 1.0,
                    'tempfactor': float(line[60:66]) if len(line) > 60 else 0.0,
                    'element': line[76:78].strip() if len(line) > 76 else ''
                }
                atoms.append(atom)
    return atoms


def write_pdb(atoms, output_path, remarks=None):
    """Write atoms to PDB file"""
    with open(output_path, 'w') as f:
        if remarks:
            for remark in remarks:
                f.write(f"REMARK  {remark}\n")
        
        for i, atom in enumerate(atoms):
            serial = i + 1
            record = atom.get('record', 'ATOM')
            name = atom['name']
            resname = atom['resname']
            chain = atom.get('chain', 'A')
            resseq = atom['resseq']
            x, y, z = atom['x'], atom['y'], atom['z']
            occupancy = atom.get('occupancy', 1.0)
            tempfactor = atom.get('tempfactor', 0.0)
            element = atom.get('element', name[0])
            
            # Format atom name (left-justify if 4 chars, otherwise center)
            if len(name) == 4:
                name_fmt = name
            else:
                name_fmt = f" {name:<3s}"
            
            line = f"{record:<6s}{serial:5d} {name_fmt:4s} {resname:3s} {chain:1s}{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{tempfactor:6.2f}          {element:>2s}\n"
            f.write(line)
        
        f.write("END\n")


def load_glycan_template():
    """Load Man5 glycan template from 1GYA"""
    print("\n📦 Loading glycan template...")
    atoms = parse_pdb(GLYCAN_TEMPLATE)
    
    # Get unique residues
    residues = set((a['resname'], a['resseq']) for a in atoms)
    print(f"   Atoms: {len(atoms)}")
    print(f"   Residues: {sorted(residues)}")
    
    # Find C1 of first NAG (attachment point)
    c1_atom = None
    for atom in atoms:
        if atom['resname'] == 'NAG' and atom['resseq'] == 1 and atom['name'] == 'C1':
            c1_atom = atom
            break
    
    if c1_atom:
        print(f"   NAG-1 C1 position: ({c1_atom['x']:.2f}, {c1_atom['y']:.2f}, {c1_atom['z']:.2f})")
    
    return atoms, c1_atom


def load_glut1():
    """Load GLUT1 and find residue 45"""
    print("\n🔬 Loading GLUT1 (4PYP)...")
    
    atoms = parse_pdb(RAW_PDB)
    
    # Filter to protein only (remove water, ligands)
    protein_atoms = [a for a in atoms if a['record'] == 'ATOM']
    print(f"   Protein atoms: {len(protein_atoms)}")
    
    # Find residue 45 (THR in crystal, should be ASN)
    res45_atoms = [a for a in protein_atoms if a['resseq'] == 45]
    print(f"   Residue 45 atoms: {len(res45_atoms)}")
    
    if res45_atoms:
        res_name = res45_atoms[0]['resname']
        print(f"   Residue 45 type: {res_name}")
        
        # Find CA and CB positions
        ca_atom = next((a for a in res45_atoms if a['name'] == 'CA'), None)
        cb_atom = next((a for a in res45_atoms if a['name'] == 'CB'), None)
        
        if ca_atom:
            print(f"   CA position: ({ca_atom['x']:.2f}, {ca_atom['y']:.2f}, {ca_atom['z']:.2f})")
    
    return protein_atoms, res45_atoms


def calculate_outward_direction(protein_atoms, res45_atoms):
    """Calculate direction pointing outward from protein surface at res45"""
    # Get protein center of mass
    coords = np.array([[a['x'], a['y'], a['z']] for a in protein_atoms])
    protein_com = np.mean(coords, axis=0)
    
    # Get res45 CA position
    ca_atom = next((a for a in res45_atoms if a['name'] == 'CA'), None)
    if ca_atom is None:
        return np.array([0, 1, 0])  # Default upward
    
    res45_pos = np.array([ca_atom['x'], ca_atom['y'], ca_atom['z']])
    
    # Direction from protein center to res45 (outward)
    outward = res45_pos - protein_com
    outward = outward / np.linalg.norm(outward)
    
    print(f"\n📐 Outward direction: ({outward[0]:.3f}, {outward[1]:.3f}, {outward[2]:.3f})")
    
    return outward


def transform_glycan(glycan_atoms, c1_atom, target_pos, outward_dir):
    """
    Transform glycan to attach at target position, pointing outward.
    
    1. Translate so C1 is at origin
    2. Rotate to align glycan axis with outward direction
    3. Translate to target position
    """
    print("\n🔄 Transforming glycan...")
    
    # Get glycan coordinates
    coords = np.array([[a['x'], a['y'], a['z']] for a in glycan_atoms])
    
    # Step 1: Translate C1 to origin
    c1_pos = np.array([c1_atom['x'], c1_atom['y'], c1_atom['z']])
    coords = coords - c1_pos
    
    # Step 2: Calculate glycan principal axis (from C1 toward glycan center)
    glycan_com = np.mean(coords, axis=0)
    glycan_axis = glycan_com / np.linalg.norm(glycan_com) if np.linalg.norm(glycan_com) > 0.1 else np.array([0, 0, 1])
    
    # Calculate rotation to align glycan_axis with outward_dir
    # Using Rodrigues' rotation formula
    v = np.cross(glycan_axis, outward_dir)
    c = np.dot(glycan_axis, outward_dir)
    
    if np.linalg.norm(v) > 1e-6:
        s = np.linalg.norm(v)
        v = v / s
        
        # Rotation matrix
        vx = np.array([[0, -v[2], v[1]], 
                       [v[2], 0, -v[0]], 
                       [-v[1], v[0], 0]])
        R = np.eye(3) + vx * s + vx @ vx * (1 - c)
        
        coords = coords @ R.T
    
    # Step 3: Translate to target position (slightly offset along outward direction)
    # Place C1 about 1.5Å from target (typical N-C bond length)
    attachment_offset = outward_dir * 1.5
    coords = coords + target_pos + attachment_offset
    
    # Update atom coordinates
    transformed_atoms = []
    for i, atom in enumerate(glycan_atoms):
        new_atom = atom.copy()
        new_atom['x'] = coords[i, 0]
        new_atom['y'] = coords[i, 1]
        new_atom['z'] = coords[i, 2]
        new_atom['chain'] = 'A'  # Same chain as protein
        # Renumber residues to avoid conflicts
        new_atom['resseq'] = atom['resseq'] + 500  # Offset residue numbers
        transformed_atoms.append(new_atom)
    
    # Report new C1 position
    new_c1 = coords[0]
    print(f"   New C1 position: ({new_c1[0]:.2f}, {new_c1[1]:.2f}, {new_c1[2]:.2f})")
    
    return transformed_atoms


def mutate_thr45_to_asn(protein_atoms):
    """
    Mutate THR45 to ASN45 (restore original glycosylation site).
    THR has: N, CA, C, O, CB, OG1, CG2
    ASN has: N, CA, C, O, CB, CG, OD1, ND2
    
    We'll keep backbone and CB, then add CG, OD1, ND2 based on CB position.
    """
    print("\n🧬 Mutating THR45 → ASN45...")
    
    new_atoms = []
    res45_backbone = []
    cb_pos = None
    ca_pos = None
    
    for atom in protein_atoms:
        if atom['resseq'] == 45:
            # Keep backbone atoms (N, CA, C, O) and CB
            if atom['name'] in ['N', 'CA', 'C', 'O', 'CB']:
                new_atom = atom.copy()
                new_atom['resname'] = 'ASN'
                new_atoms.append(new_atom)
                
                if atom['name'] == 'CB':
                    cb_pos = np.array([atom['x'], atom['y'], atom['z']])
                if atom['name'] == 'CA':
                    ca_pos = np.array([atom['x'], atom['y'], atom['z']])
            # Skip THR-specific atoms (OG1, CG2, HG1, etc.)
        else:
            new_atoms.append(atom)
    
    # Add ASN sidechain atoms (CG, OD1, ND2)
    if cb_pos is not None and ca_pos is not None:
        # Direction from CA to CB
        cb_dir = cb_pos - ca_pos
        cb_dir = cb_dir / np.linalg.norm(cb_dir)
        
        # CG is ~1.5Å from CB along similar direction
        cg_pos = cb_pos + cb_dir * 1.52
        
        # OD1 and ND2 are ~1.2-1.3Å from CG, roughly perpendicular
        perp = np.cross(cb_dir, np.array([0, 0, 1]))
        if np.linalg.norm(perp) < 0.1:
            perp = np.cross(cb_dir, np.array([0, 1, 0]))
        perp = perp / np.linalg.norm(perp)
        
        od1_pos = cg_pos + perp * 1.23
        nd2_pos = cg_pos - perp * 1.32  # ND2 is where glycan attaches
        
        # Find insertion point (after CB of res45)
        insert_idx = None
        for i, atom in enumerate(new_atoms):
            if atom['resseq'] == 45 and atom['name'] == 'CB':
                insert_idx = i + 1
                break
        
        if insert_idx:
            asn_atoms = [
                {'record': 'ATOM', 'name': 'CG', 'resname': 'ASN', 'chain': 'A', 
                 'resseq': 45, 'x': cg_pos[0], 'y': cg_pos[1], 'z': cg_pos[2],
                 'occupancy': 1.0, 'tempfactor': 0.0, 'element': 'C'},
                {'record': 'ATOM', 'name': 'OD1', 'resname': 'ASN', 'chain': 'A',
                 'resseq': 45, 'x': od1_pos[0], 'y': od1_pos[1], 'z': od1_pos[2],
                 'occupancy': 1.0, 'tempfactor': 0.0, 'element': 'O'},
                {'record': 'ATOM', 'name': 'ND2', 'resname': 'ASN', 'chain': 'A',
                 'resseq': 45, 'x': nd2_pos[0], 'y': nd2_pos[1], 'z': nd2_pos[2],
                 'occupancy': 1.0, 'tempfactor': 0.0, 'element': 'N'},
            ]
            
            for i, asn_atom in enumerate(asn_atoms):
                new_atoms.insert(insert_idx + i, asn_atom)
            
            print(f"   Added ASN sidechain: CG, OD1, ND2")
            print(f"   ND2 position (glycan attachment): ({nd2_pos[0]:.2f}, {nd2_pos[1]:.2f}, {nd2_pos[2]:.2f})")
            
            return new_atoms, nd2_pos
    
    return new_atoms, None


def main():
    print("=" * 70)
    print("Phase 4: Building Glycosylated GLUT1 (Template-based)")
    print("=" * 70)
    print("\nUsing Man5 N-glycan template from PDB 1GYA")
    
    # Step 1: Load glycan template
    glycan_atoms, c1_atom = load_glycan_template()
    
    if c1_atom is None:
        print("❌ Could not find C1 atom in glycan template")
        return
    
    # Step 2: Load GLUT1
    protein_atoms, res45_atoms = load_glut1()
    
    # Step 3: Mutate THR45 → ASN45
    protein_atoms, nd2_pos = mutate_thr45_to_asn(protein_atoms)
    
    if nd2_pos is None:
        print("❌ Could not create ASN45")
        return
    
    # Step 4: Calculate outward direction
    outward_dir = calculate_outward_direction(protein_atoms, 
                                               [a for a in protein_atoms if a['resseq'] == 45])
    
    # Step 5: Transform glycan to attachment site
    transformed_glycan = transform_glycan(glycan_atoms, c1_atom, nd2_pos, outward_dir)
    
    # Step 6: Save naked GLUT1 (with ASN45 mutation)
    naked_pdb = os.path.join(STRUCTURES_DIR, "glut1_naked_asn45.pdb")
    write_pdb(protein_atoms, naked_pdb, 
              remarks=["GLUT1 with THR45 mutated to ASN45 (glycosylation site)",
                       "No glycan attached (naked model for cancer cell)"])
    print(f"\n✅ Naked GLUT1 saved: {naked_pdb}")
    
    # Step 7: Merge protein + glycan
    merged_atoms = protein_atoms + transformed_glycan
    
    glyco_pdb = os.path.join(STRUCTURES_DIR, "glut1_glycosylated_man5.pdb")
    write_pdb(merged_atoms, glyco_pdb,
              remarks=["Glycosylated GLUT1 - Man5 N-glycan at ASN45",
                       "Glycan template from PDB 1GYA",
                       "Normal cell model (RBC, Endothelial)"])
    print(f"✅ Glycosylated GLUT1 saved: {glyco_pdb}")
    
    # Summary
    print("\n" + "=" * 70)
    print("GLYCAN STRUCTURE SUMMARY")
    print("=" * 70)
    
    glycan_residues = set((a['resname'], a['resseq']) for a in transformed_glycan)
    print(f"\nGlycan composition:")
    for resname, resseq in sorted(glycan_residues):
        count = len([a for a in transformed_glycan if a['resseq'] == resseq])
        print(f"   {resname} {resseq}: {count} atoms")
    
    print(f"\nTotal glycan atoms: {len(transformed_glycan)}")
    print(f"Total protein atoms: {len(protein_atoms)}")
    print(f"Total merged atoms: {len(merged_atoms)}")
    
    # Calculate glycan extent
    glycan_coords = np.array([[a['x'], a['y'], a['z']] for a in transformed_glycan])
    glycan_extent = np.max(glycan_coords, axis=0) - np.min(glycan_coords, axis=0)
    print(f"\nGlycan dimensions: {glycan_extent[0]:.1f} x {glycan_extent[1]:.1f} x {glycan_extent[2]:.1f} Å")
    
    print("\n" + "=" * 70)
    print("Model Building Complete!")
    print("=" * 70)
    print(f"\nOutput files:")
    print(f"   1. Naked GLUT1: {naked_pdb}")
    print(f"   2. Glycosylated GLUT1: {glyco_pdb}")


if __name__ == "__main__":
    main()
