#!/usr/bin/env python
"""
DiffDock Docking for Short-Tripod Models
=========================================
Uses deep learning-based diffusion model for blind docking
"""

import os
import sys
import subprocess
import shutil
from rdkit import Chem

# Paths
DIFFDOCK_DIR = "/home/pjho3/projects/Drug/DiffDock"
RECEPTOR_PDB = "/home/pjho3/projects/Drug/raw_data/4PYP.pdb"
STRUCTURES_DIR = "/home/pjho3/projects/Drug/structures"
OUTPUT_DIR = "/home/pjho3/projects/Drug/results/diffdock_docking"

os.makedirs(OUTPUT_DIR, exist_ok=True)

def get_smiles_from_sdf(sdf_path):
    """Extract SMILES from SDF file"""
    mol = Chem.SDMolSupplier(sdf_path, removeHs=False)[0]
    if mol is None:
        raise ValueError(f"Failed to load {sdf_path}")
    return Chem.MolToSmiles(Chem.RemoveHs(mol))

def run_diffdock(receptor_pdb, ligand_smiles, output_dir, complex_name):
    """Run DiffDock inference"""
    print(f"\n{'='*60}")
    print(f"DiffDock Docking: {complex_name}")
    print(f"{'='*60}")
    
    # Create output directory for this complex
    complex_out = os.path.join(output_dir, complex_name)
    os.makedirs(complex_out, exist_ok=True)
    
    # DiffDock command
    cmd = [
        "mamba", "run", "-n", "diffdock",
        "python", "inference.py",
        "--protein_path", receptor_pdb,
        "--ligand_description", ligand_smiles,
        "--out_dir", complex_out,
        "--complex_name", complex_name,
        "--samples_per_complex", "10",
        "--inference_steps", "20",
        "--batch_size", "10"
    ]
    
    print(f"SMILES: {ligand_smiles[:50]}...")
    print("Running DiffDock (this may take a few minutes)...")
    
    result = subprocess.run(
        cmd,
        cwd=DIFFDOCK_DIR,
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        print(f"❌ DiffDock failed!")
        print(f"STDERR: {result.stderr[-1000:]}")
        return None
    
    # Find output files
    output_files = []
    for f in os.listdir(complex_out):
        if f.endswith('.sdf'):
            output_files.append(os.path.join(complex_out, f))
    
    if output_files:
        print(f"✅ Generated {len(output_files)} pose files")
        
        # Parse confidence scores from filenames
        print("\nTop poses:")
        scored_files = []
        for f in output_files:
            basename = os.path.basename(f)
            if 'confidence' in basename:
                try:
                    conf = float(basename.split('confidence')[-1].replace('.sdf', ''))
                    scored_files.append((conf, f))
                except:
                    pass
        
        scored_files.sort(reverse=True)
        for i, (conf, path) in enumerate(scored_files[:5]):
            print(f"  Pose {i+1}: confidence={conf:.3f}")
        
        return complex_out
    else:
        print(f"❌ No output generated")
        return None

def main():
    print("="*60)
    print("DiffDock Docking for Short-Tripod Models")
    print("="*60)
    
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
        
        # Get SMILES
        try:
            smiles = get_smiles_from_sdf(ligand_path)
        except Exception as e:
            print(f"❌ Failed to get SMILES: {e}")
            continue
        
        output = run_diffdock(
            RECEPTOR_PDB,
            smiles,
            OUTPUT_DIR,
            model_name
        )
        
        if output:
            results[model_name] = output
    
    print("\n" + "="*60)
    print("DiffDock Docking Complete!")
    print("="*60)
    print(f"Output directory: {OUTPUT_DIR}")
    for name, path in results.items():
        print(f"  {name}: {path}")

if __name__ == "__main__":
    main()
