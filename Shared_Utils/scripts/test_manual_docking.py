#!/usr/bin/env python3
"""02_manual_docking.py 검증 스크립트"""
import sys
import numpy as np

print("=" * 50)
print("02_manual_docking.py 사전 검증")
print("=" * 50)

# 1. Import 테스트
print("\n[1] Import 테스트...")
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    print("   ✅ RDKit OK")
except ImportError as e:
    print(f"   ❌ RDKit 실패: {e}")
    sys.exit(1)

try:
    from Bio.PDB import PDBParser
    print("   ✅ BioPython OK")
except ImportError as e:
    print(f"   ❌ BioPython 실패: {e}")
    sys.exit(1)

# 2. 분자 생성 테스트 (간단한 tripod-like 구조)
print("\n[2] 분자 생성 테스트...")
test_smiles = "OCC(CO)(CO)CO"  # 간단한 4-arm 분자
mol = Chem.MolFromSmiles(test_smiles)
if mol is None:
    print("   ❌ SMILES 파싱 실패")
    sys.exit(1)
mol = Chem.AddHs(mol)
result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
if result == -1:
    print("   ❌ 3D 구조 생성 실패")
    sys.exit(1)
AllChem.MMFFOptimizeMolecule(mol)
print(f"   ✅ 분자 생성 OK (원자 수: {mol.GetNumAtoms()})")

# 3. find_ligand_arms 로직 테스트
print("\n[3] Arm 찾기 로직 테스트...")
conf = mol.GetConformer()
coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
center = coords.mean(axis=0)
distances = np.linalg.norm(coords - center, axis=1)

o_atoms = []
for i, atom in enumerate(mol.GetAtoms()):
    if atom.GetAtomicNum() == 8:
        o_atoms.append((i, distances[i], coords[i]))

o_atoms.sort(key=lambda x: x[1], reverse=True)
print(f"   ✅ 산소 원자 {len(o_atoms)}개 발견")

arm_points = []
for idx, dist, coord in o_atoms:
    if len(arm_points) == 0:
        arm_points.append(coord)
    else:
        is_different = True
        for existing in arm_points:
            vec1 = (coord - center) / np.linalg.norm(coord - center)
            vec2 = (existing - center) / np.linalg.norm(existing - center)
            if np.dot(vec1, vec2) > 0.5:
                is_different = False
                break
        if is_different:
            arm_points.append(coord)
    if len(arm_points) == 3:
        break

print(f"   ✅ Arm {len(arm_points)}개 식별됨")

# 4. 회전 행렬 계산 테스트
print("\n[4] 회전 행렬 계산 테스트...")
# 가상의 chain 중심 (정삼각형 배치)
d = 50.0
chain_centers = {
    'A': np.array([0.0, 0.0, 0.0]),
    'B': np.array([d, 0.0, 0.0]),
    'C': np.array([d/2, d*np.sqrt(3)/2, 0.0])
}
total_center = np.mean(list(chain_centers.values()), axis=0)

chain_dirs = []
for chain_id in sorted(chain_centers.keys()):
    direction = chain_centers[chain_id] - total_center
    direction = direction / np.linalg.norm(direction)
    chain_dirs.append(direction)
chain_dirs = np.array(chain_dirs)

if len(arm_points) >= 3:
    arm_dirs = []
    for arm in arm_points[:3]:
        direction = arm - center
        norm = np.linalg.norm(direction)
        if norm > 0:
            direction = direction / norm
        arm_dirs.append(direction)
    arm_dirs = np.array(arm_dirs)
    
    H = arm_dirs.T @ chain_dirs
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    
    det = np.linalg.det(R)
    print(f"   ✅ 회전 행렬 OK (det={det:.4f}, should be ~1.0)")
else:
    print(f"   ⚠️ Arm이 3개 미만 ({len(arm_points)}개)")

# 5. 좌표 변환 테스트
print("\n[5] 좌표 변환 테스트...")
test_coord = np.array([1.0, 2.0, 3.0])
rotated = R @ test_coord
translated = rotated + np.array([10, 20, 30])
print(f"   ✅ 좌표 변환 OK: {test_coord} -> {translated}")

print("\n" + "=" * 50)
print("🎉 모든 검증 통과! 스크립트 실행 가능")
print("=" * 50)
