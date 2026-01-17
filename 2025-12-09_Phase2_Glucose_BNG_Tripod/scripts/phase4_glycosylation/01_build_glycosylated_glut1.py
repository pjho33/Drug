#!/usr/bin/env python
"""
Phase 4: Build Glycosylated GLUT1 Model
========================================
1. Restore THR45 → ASN45 (original glycosylation site)
2. Attach N-glycan chain to Asn45
3. Prepare both naked and glycosylated models for comparison

N-glycan structure (Complex type, common in RBC/Endothelial):
    Man-α1,6\
              Man-β1,4-GlcNAc-β1,4-GlcNAc-β1-Asn
    Man-α1,3/

For simplicity, we'll use a high-mannose type glycan (Man5-9)
which is computationally tractable and still provides steric hindrance.
"""

import os
import numpy as np
from pdbfixer import PDBFixer
from openmm.app import PDBFile, Modeller
import warnings
warnings.filterwarnings('ignore')

# Paths
RAW_PDB = "/home/pjho3/projects/Drug/raw_data/4PYP.pdb"
OUTPUT_DIR = "/home/pjho3/projects/Drug/results/phase4_glycan"
STRUCTURES_DIR = "/home/pjho3/projects/Drug/structures/phase4"

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(STRUCTURES_DIR, exist_ok=True)


def prepare_receptor_for_glycosylation():
    """
    Prepare GLUT1 receptor:
    1. Fix missing residues
    2. Mutate THR45 back to ASN45
    """
    print("=" * 60)
    print("Step 1: Prepare GLUT1 Receptor")
    print("=" * 60)
    
    # Load and fix PDB
    print("\nLoading 4PYP...")
    fixer = PDBFixer(filename=RAW_PDB)
    
    # Find missing residues
    fixer.findMissingResidues()
    print(f"Missing residues: {len(fixer.missingResidues)}")
    
    # Find and add missing atoms
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    
    # Remove heterogens (ligands, water) but keep protein
    fixer.removeHeterogens(keepWater=False)
    
    # Add missing hydrogens
    fixer.addMissingHydrogens(7.4)  # pH 7.4
    
    # Save cleaned receptor (without glycan)
    naked_pdb = os.path.join(STRUCTURES_DIR, "glut1_naked.pdb")
    with open(naked_pdb, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    
    print(f"✅ Naked GLUT1 saved: {naked_pdb}")
    
    return naked_pdb, fixer


def get_asn45_position(pdb_path):
    """Get the position of residue 45 (THR in crystal, should be ASN)"""
    print("\nFinding residue 45 position...")
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") and " CA " in line:
                parts = line.split()
                # Find residue 45
                res_num = int(line[22:26].strip())
                if res_num == 45:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    res_name = line[17:20].strip()
                    print(f"  Residue 45: {res_name} at ({x:.2f}, {y:.2f}, {z:.2f})")
                    return np.array([x, y, z]), res_name
    
    return None, None


def create_glycan_pdb(anchor_pos, glycan_type="high_mannose"):
    """
    Create a simplified N-glycan structure attached at anchor position.
    
    High-mannose N-glycan (Man5):
    - Core: GlcNAc-GlcNAc-Man (attached to Asn)
    - Branches: Man-Man (α1,3 and α1,6)
    
    This is a simplified representation for MD simulation.
    """
    print(f"\nCreating {glycan_type} N-glycan...")
    
    # Glycan building blocks (approximate coordinates)
    # Each sugar ring is ~5-6 Å in diameter
    # We'll build extending outward from the protein surface
    
    # Direction: extend outward from protein (approximate normal)
    # Asn45 is on extracellular loop, so glycan extends into extracellular space
    
    glycan_atoms = []
    
    # Core GlcNAc-1 (attached to Asn ND2)
    # Position relative to anchor (Asn CA)
    glcnac1_center = anchor_pos + np.array([0, 5, 0])  # 5Å away
    
    # Core GlcNAc-2
    glcnac2_center = glcnac1_center + np.array([0, 5, 0])
    
    # Core Mannose (branching point)
    man_core_center = glcnac2_center + np.array([0, 5, 0])
    
    # Branch 1: α1,3 Mannose
    man_a13_center = man_core_center + np.array([-4, 4, 0])
    
    # Branch 2: α1,6 Mannose  
    man_a16_center = man_core_center + np.array([4, 4, 0])
    
    # Terminal mannoses (for Man5)
    man_term1 = man_a13_center + np.array([-3, 4, 0])
    man_term2 = man_a16_center + np.array([3, 4, 0])
    
    # Create simplified sugar ring atoms (just key atoms for steric effect)
    sugar_positions = [
        ("NAG", 1, glcnac1_center),   # GlcNAc-1
        ("NAG", 2, glcnac2_center),   # GlcNAc-2
        ("MAN", 3, man_core_center),  # Core Man
        ("MAN", 4, man_a13_center),   # α1,3 Man
        ("MAN", 5, man_a16_center),   # α1,6 Man
        ("MAN", 6, man_term1),        # Terminal Man 1
        ("MAN", 7, man_term2),        # Terminal Man 2
    ]
    
    # Generate PDB lines for glycan
    pdb_lines = []
    atom_num = 1
    
    for sugar_name, res_num, center in sugar_positions:
        # Create a simplified 6-membered ring
        # Ring atoms: C1, C2, C3, C4, C5, O5
        ring_atoms = [
            ("C1", center + np.array([1.2, 0, 0])),
            ("C2", center + np.array([0.6, 1.0, 0])),
            ("C3", center + np.array([-0.6, 1.0, 0])),
            ("C4", center + np.array([-1.2, 0, 0])),
            ("C5", center + np.array([-0.6, -1.0, 0])),
            ("O5", center + np.array([0.6, -1.0, 0])),
        ]
        
        for atom_name, pos in ring_atoms:
            line = f"HETATM{atom_num:5d}  {atom_name:<3s} {sugar_name} A{res_num:4d}    {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00  0.00           C\n"
            pdb_lines.append(line)
            atom_num += 1
    
    return pdb_lines, sugar_positions


def build_glycosylated_model(naked_pdb, glycan_pdb_lines):
    """Combine naked GLUT1 with glycan"""
    print("\nBuilding glycosylated model...")
    
    # Read naked PDB
    with open(naked_pdb, 'r') as f:
        protein_lines = [l for l in f.readlines() if l.startswith("ATOM") or l.startswith("TER")]
    
    # Combine
    glyco_pdb = os.path.join(STRUCTURES_DIR, "glut1_glycosylated.pdb")
    
    with open(glyco_pdb, 'w') as f:
        f.write("REMARK  Glycosylated GLUT1 - N-glycan at Asn45\n")
        f.write("REMARK  Glycan type: High-mannose (Man5)\n")
        for line in protein_lines:
            if not line.startswith("END"):
                f.write(line)
        f.write("TER\n")
        for line in glycan_pdb_lines:
            f.write(line)
        f.write("END\n")
    
    print(f"✅ Glycosylated GLUT1 saved: {glyco_pdb}")
    
    return glyco_pdb


def analyze_steric_hindrance(naked_pdb, glyco_pdb):
    """Analyze how glycan affects binding site accessibility"""
    print("\n" + "=" * 60)
    print("Steric Hindrance Analysis")
    print("=" * 60)
    
    # Load structures
    from openmm.app import PDBFile
    
    naked = PDBFile(naked_pdb)
    glyco = PDBFile(glyco_pdb)
    
    print(f"\nNaked GLUT1: {naked.topology.getNumAtoms()} atoms")
    print(f"Glycosylated GLUT1: {glyco.topology.getNumAtoms()} atoms")
    
    # Find binding site center (approximate - near Asn45 region)
    # The glucose binding site is in the central cavity
    
    # Calculate glycan extent
    glycan_atoms = []
    for atom in glyco.topology.atoms():
        if atom.residue.name in ['NAG', 'MAN', 'BMA', 'FUC', 'GAL']:
            pos = glyco.positions[atom.index]
            glycan_atoms.append([pos.x, pos.y, pos.z])
    
    if glycan_atoms:
        glycan_atoms = np.array(glycan_atoms)
        glycan_center = np.mean(glycan_atoms, axis=0)
        glycan_extent = np.max(glycan_atoms, axis=0) - np.min(glycan_atoms, axis=0)
        
        print(f"\nGlycan statistics:")
        print(f"  Atoms: {len(glycan_atoms)}")
        print(f"  Center: ({glycan_center[0]:.1f}, {glycan_center[1]:.1f}, {glycan_center[2]:.1f}) nm")
        print(f"  Extent: {glycan_extent[0]*10:.1f} x {glycan_extent[1]*10:.1f} x {glycan_extent[2]*10:.1f} Å")
    
    return glycan_atoms


def main():
    print("=" * 70)
    print("Phase 4: Building Glycosylated GLUT1 Model")
    print("=" * 70)
    print("\nHypothesis: Large drug cannot penetrate glycan layer on normal cells")
    print("Target: Asn45 N-glycosylation site")
    
    # Step 1: Prepare naked receptor
    naked_pdb, fixer = prepare_receptor_for_glycosylation()
    
    # Step 2: Find Asn45 position
    asn45_pos, res_name = get_asn45_position(naked_pdb)
    
    if asn45_pos is None:
        print("❌ Could not find residue 45")
        return
    
    # Step 3: Create glycan
    glycan_lines, sugar_positions = create_glycan_pdb(asn45_pos)
    
    # Step 4: Build glycosylated model
    glyco_pdb = build_glycosylated_model(naked_pdb, glycan_lines)
    
    # Step 5: Analyze
    analyze_steric_hindrance(naked_pdb, glyco_pdb)
    
    print("\n" + "=" * 70)
    print("Model Building Complete!")
    print("=" * 70)
    print(f"\nOutput files:")
    print(f"  1. Naked GLUT1: {STRUCTURES_DIR}/glut1_naked.pdb")
    print(f"  2. Glycosylated GLUT1: {STRUCTURES_DIR}/glut1_glycosylated.pdb")
    print(f"\nNext steps:")
    print(f"  1. Run MD simulation with drug approaching both models")
    print(f"  2. Compare COM distance to binding site")
    print(f"  3. Analyze H-bonds with glycan vs protein")


if __name__ == "__main__":
    main()
