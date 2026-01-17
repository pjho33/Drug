# scripts/09_analyze_size.py
"""Analyze molecular size of Tripod and control ligands."""
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import numpy as np

# SMILES definitions
MOLECULES = {
    "Tripod (PEG2-L-glucose x3)": "NC(COCCOCCOOC1OC(CO)C(O)C(O)C(O)C1O)(COCCOCCOOC1OC(CO)C(O)C(O)C(O)C1O)COCCOCCOOC1OC(CO)C(O)C(O)C(O)C1O",
    "D-glucose": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",
    "BNG": "CCCCCCCCCO[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"
}

def analyze_molecule(name, smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates with multiple attempts
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.maxAttempts = 1000
    result = AllChem.EmbedMolecule(mol, params)
    
    if result == -1:
        # Try with more permissive parameters
        params2 = AllChem.ETKDG()
        params2.randomSeed = 42
        params2.useRandomCoords = True
        params2.maxAttempts = 5000
        result = AllChem.EmbedMolecule(mol, params2)
    
    if result != -1:
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        except:
            try:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
            except:
                pass
    
    # Get conformer
    conf = mol.GetConformer()
    positions = conf.GetPositions()
    
    # Calculate molecular properties
    mw = Descriptors.MolWt(mol)
    n_heavy = mol.GetNumHeavyAtoms()
    
    # Calculate size metrics
    # 1. Maximum distance between any two atoms (diameter)
    max_dist = 0
    for i in range(len(positions)):
        for j in range(i+1, len(positions)):
            dist = np.linalg.norm(positions[i] - positions[j])
            if dist > max_dist:
                max_dist = dist
    
    # 2. Radius of gyration
    center = np.mean(positions, axis=0)
    rg = np.sqrt(np.mean(np.sum((positions - center)**2, axis=1)))
    
    # 3. Bounding box dimensions
    min_coords = np.min(positions, axis=0)
    max_coords = np.max(positions, axis=0)
    dimensions = max_coords - min_coords
    volume = dimensions[0] * dimensions[1] * dimensions[2]
    
    return {
        "name": name,
        "mw": mw,
        "n_heavy": n_heavy,
        "diameter": max_dist,
        "rg": rg,
        "dimensions": dimensions,
        "volume": volume
    }

if __name__ == "__main__":
    print("\n" + "#"*60)
    print("MOLECULAR SIZE COMPARISON")
    print("#"*60)
    
    results = []
    for name, smiles in MOLECULES.items():
        data = analyze_molecule(name, smiles)
        results.append(data)
        
        print(f"\n{'='*50}")
        print(f"{name}")
        print(f"{'='*50}")
        print(f"Molecular Weight: {data['mw']:.1f} Da")
        print(f"Heavy Atoms: {data['n_heavy']}")
        print(f"\nSize Metrics:")
        print(f"  Max Diameter: {data['diameter']:.1f} A")
        print(f"  Radius of Gyration: {data['rg']:.1f} A")
        dim = data['dimensions']
        print(f"  Bounding Box: {dim[0]:.1f} x {dim[1]:.1f} x {dim[2]:.1f} A")
        print(f"  Approx Volume: {data['volume']:.0f} A^3")
    
    # Summary table
    print("\n" + "="*70)
    print("SUMMARY COMPARISON")
    print("="*70)
    print(f"{'Molecule':<30} {'MW (Da)':<12} {'Diameter (A)':<15} {'Rg (A)':<10}")
    print("-"*70)
    for r in results:
        print(f"{r['name']:<30} {r['mw']:<12.1f} {r['diameter']:<15.1f} {r['rg']:<10.1f}")
    print("="*70)
    
    # Comparison ratios
    tripod = results[0]
    glucose = results[1]
    print(f"\nTripod vs D-glucose:")
    print(f"  - {tripod['mw']/glucose['mw']:.1f}x heavier")
    print(f"  - {tripod['diameter']/glucose['diameter']:.1f}x larger diameter")
    print(f"  - {tripod['volume']/glucose['volume']:.1f}x larger volume")
