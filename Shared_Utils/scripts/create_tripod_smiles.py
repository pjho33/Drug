# scripts/create_tripod_smiles.py
# Tripod 리간드 SMILES 생성 - Center에 triazole core, PEG arm 끝에 L-glucose

from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

print("=" * 60)
print("Tripod SMILES Generator")
print("Structure: Benzene-CH2-Triazole-CH2-PEG6-L-glucose")
print("=" * 60)

# 올바른 구조:
# Center: 1,3,5-trisubstituted benzene (rigid core)
# Linker: -CH2-[1,2,3-triazole]-CH2-PEG6-
# End: L-glucose (at the END of PEG, not at center!)

# Step 1: L-glucose 확인
# L-glucose (alpha-L-glucopyranose) - C1에서 O-glycosidic linkage
l_glucose = "O[C@H]1[C@H](O)[C@@H](O)[C@H](CO)O[C@H]1O"
mol_glc = Chem.MolFromSmiles(l_glucose)
print(f"L-glucose atoms: {mol_glc.GetNumAtoms() if mol_glc else 'Invalid'}")

# Step 2: PEG6 linker (6 ethylene glycol units)
# Structure: -CH2-CH2-O-CH2-CH2-O-CH2-CH2-O-CH2-CH2-O-CH2-CH2-O-CH2-CH2-O-
peg6 = "CCOCCOCCOCCOCCOCCO"

# Step 3: 1,2,3-triazole (1,4-disubstituted)
# N1 position: connected to benzyl (Bn-CH2-N)
# C4 position: connected to PEG linker
# SMILES: n1cc(R)nn1 where R is at C4

# Step 4: Build the arm
# Benzene-CH2-N(triazole C4)-CH2-PEG6-O-L-glucose
# The triazole: Cn1ccnn1 (N1-substituted), need C4 substitution

# Correct 1,4-disubstituted 1,2,3-triazole:
# n1cc(substituent)nn1 - this puts substituent at C4
# Cn1cc(...)nn1 - N1 has CH2 from benzene, C4 has the PEG chain

# L-glucose connected via O-glycosidic bond at C1 (anomeric)
# PEG-O-C1(glucose)
# Final glucose part: O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O

# Build arm: triazole-CH2-PEG6-O-glucose
# Cn1cc(COCCOCCOCCOCCOCCO[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)nn1

# L-glucose with glycosidic linkage (beta form for stability)
l_glc_linked = "O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"

# Full arm
arm = f"Cn1cc(COCCOCCOCCOCCOCCO{l_glc_linked})nn1"

# Test arm first
mol_arm = Chem.MolFromSmiles(arm)
if mol_arm:
    print(f"✅ Single arm valid, atoms: {mol_arm.GetNumAtoms()}")
else:
    print("❌ Arm invalid, trying alternative...")
    # Alternative L-glucose linkage
    l_glc_linked = "[C@@H]1(O)O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"
    arm = f"Cn1cc(COCCOCCOCCOCCOCCO{l_glc_linked})nn1"
    mol_arm = Chem.MolFromSmiles(arm)
    print(f"Alternative arm: {'valid' if mol_arm else 'invalid'}")

# Build full tripod: 1,3,5-trisubstituted benzene
tripod_smiles = f"c1c({arm})cc({arm})cc1{arm}"

mol = Chem.MolFromSmiles(tripod_smiles)

if mol:
    print(f"\n✅ Full Tripod SMILES is valid!")
    print(f"Total atoms: {mol.GetNumAtoms()}")
    
    ring_info = mol.GetRingInfo()
    print(f"Number of rings: {ring_info.NumRings()} (expected: 7 = 1 benzene + 3 triazole + 3 glucose)")
    
    # 파일로 저장
    with open("/home/pjho3/projects/Drug/tripod_peg6_l_glucose_v2.smi", "w") as f:
        f.write(tripod_smiles)
    print(f"\n✅ Saved to tripod_peg6_l_glucose_v2.smi")
    print(f"\nSMILES:\n{tripod_smiles}")
    
    # 3D 구조 생성
    print("\n" + "=" * 60)
    print("Generating 3D structure...")
    print("=" * 60)
    
    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.maxIterations = 10000
    params.randomSeed = 42
    params.useRandomCoords = True
    
    result = AllChem.EmbedMolecule(mol_h, params)
    if result == 0:
        print("✅ 3D embedding successful!")
        
        try:
            AllChem.MMFFOptimizeMolecule(mol_h, maxIters=2000)
            print("✅ MMFF optimization done!")
        except:
            AllChem.UFFOptimizeMolecule(mol_h, maxIters=2000)
            print("✅ UFF optimization done!")
        
        # Arm 길이 측정 (center에서 glucose까지)
        conf = mol_h.GetConformer()
        coords = np.array([conf.GetAtomPosition(i) for i in range(mol_h.GetNumAtoms())])
        center = coords.mean(axis=0)
        
        # 가장 먼 원자까지의 거리 (arm 길이 추정)
        distances = np.linalg.norm(coords - center, axis=1)
        max_dist = distances.max()
        print(f"📏 Estimated arm length (center to tip): {max_dist:.1f} Å")
        
        # SDF 저장
        writer = Chem.SDWriter("/home/pjho3/projects/Drug/tripod_peg6_l_glucose_v2.sdf")
        writer.write(mol_h)
        writer.close()
        print("✅ 3D structure saved to tripod_peg6_l_glucose_v2.sdf")
    else:
        print("❌ 3D embedding failed")
else:
    print("❌ Full Tripod SMILES is invalid")
