#!/usr/bin/env python
"""
Gnina Docking for Short-Tripod Models
======================================
Uses CNN-based scoring for accurate ligand placement
"""

import subprocess
import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

# Paths
GNINA_PATH = "/home/pjho3/projects/Drug/tools/gnina"
RECEPTOR_PDB = "/home/pjho3/projects/Drug/raw_data/4PYP.pdb"
STRUCTURES_DIR = "/home/pjho3/projects/Drug/structures"
OUTPUT_DIR = "/home/pjho3/projects/Drug/results/gnina_docking"

os.makedirs(OUTPUT_DIR, exist_ok=True)

def get_binding_site_center(pdb_path):
    """Get GLUT1 binding site center"""
    from Bio.PDB import PDBParser
    
    # GLUT1 glucose binding site residues
    BINDING_RESIDUES = [34, 161, 282, 283, 288, 317, 388, 412]
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('receptor', pdb_path)
    
    binding_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id()[1] in BINDING_RESIDUES and 'CA' in residue:
                    binding_coords.append(residue['CA'].get_coord())
    
    if binding_coords:
        center = np.array(binding_coords).mean(axis=0)
    else:
        # Fallback to geometric center
        all_coords = []
        for atom in structure.get_atoms():
            all_coords.append(atom.get_coord())
        center = np.array(all_coords).mean(axis=0)
    
    return center

def sdf_to_mol2(sdf_path, mol2_path):
    """Convert SDF to MOL2 using RDKit"""
    mol = Chem.SDMolSupplier(sdf_path, removeHs=False)[0]
    if mol is None:
        raise ValueError(f"Failed to load {sdf_path}")
    
    # Use obabel for conversion (more reliable for mol2)
    cmd = f"obabel {sdf_path} -O {mol2_path} 2>/dev/null"
    result = subprocess.run(cmd, shell=True, capture_output=True)
    
    if not os.path.exists(mol2_path):
        # Fallback: just use SDF directly with gnina
        return sdf_path
    return mol2_path

def run_gnina_docking(receptor_pdb, ligand_sdf, output_prefix, center, box_size=25):
    """Run Gnina docking"""
    print(f"\n{'='*60}")
    print(f"Gnina Docking: {os.path.basename(ligand_sdf)}")
    print(f"{'='*60}")
    
    output_sdf = f"{output_prefix}_docked.sdf"
    log_file = f"{output_prefix}_log.txt"
    
    # Gnina command
    cmd = [
        GNINA_PATH,
        "-r", receptor_pdb,
        "-l", ligand_sdf,
        "-o", output_sdf,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(box_size),
        "--size_y", str(box_size),
        "--size_z", str(box_size),
        "--exhaustiveness", "32",
        "--num_modes", "9",
        "--cnn_scoring", "rescore",  # Use CNN for rescoring
        "--seed", "42"
    ]
    
    print(f"Center: ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})")
    print(f"Box size: {box_size} Å")
    print("Running Gnina...")
    
    with open(log_file, "w") as log:
        result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    
    if result.returncode != 0:
        print(f"❌ Gnina failed! Check {log_file}")
        return None
    
    # Parse results
    if os.path.exists(output_sdf):
        suppl = Chem.SDMolSupplier(output_sdf, removeHs=False)
        mols = [m for m in suppl if m is not None]
        print(f"✅ Generated {len(mols)} poses")
        
        # Print scores
        print("\nTop poses:")
        for i, mol in enumerate(mols[:5]):
            if mol.HasProp("CNNscore"):
                cnn_score = mol.GetProp("CNNscore")
                cnn_affinity = mol.GetProp("CNNaffinity") if mol.HasProp("CNNaffinity") else "N/A"
                print(f"  Pose {i+1}: CNN score={cnn_score}, affinity={cnn_affinity}")
        
        return output_sdf
    else:
        print(f"❌ No output generated")
        return None

def main():
    print("="*60)
    print("Gnina Docking for Short-Tripod Models")
    print("="*60)
    
    # Get binding site center
    print("\nFinding binding site center...")
    center = get_binding_site_center(RECEPTOR_PDB)
    print(f"Binding site center: ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})")
    
    # Models to dock
    models = [
        ("model_a_tris_peg2", "model_a_tris_peg2.sdf"),
        ("model_b_arg_tris_peg2", "model_b_arg_tris_peg2.sdf"),
    ]
    
    results = {}
    for model_name, ligand_file in models:
        ligand_path = os.path.join(STRUCTURES_DIR, ligand_file)
        if not os.path.exists(ligand_path):
            print(f"⚠️ Ligand not found: {ligand_path}")
            continue
        
        output_prefix = os.path.join(OUTPUT_DIR, model_name)
        docked_sdf = run_gnina_docking(
            RECEPTOR_PDB, 
            ligand_path, 
            output_prefix, 
            center,
            box_size=30  # Larger box for tripod
        )
        
        if docked_sdf:
            results[model_name] = docked_sdf
    
    print("\n" + "="*60)
    print("Docking Complete!")
    print("="*60)
    print(f"Output directory: {OUTPUT_DIR}")
    for name, path in results.items():
        print(f"  {name}: {path}")

if __name__ == "__main__":
    main()
