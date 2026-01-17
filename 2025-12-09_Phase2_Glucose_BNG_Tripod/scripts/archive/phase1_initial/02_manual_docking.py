# scripts/02_manual_docking.py
import sys
import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from Bio.PDB import PDBParser

# GLUT1 glucose binding site 잔기들 (4PYP 구조 기반)
BINDING_SITE_RESIDUES = [34, 161, 282, 283, 288, 317, 388, 412]

def get_receptor_info(receptor_pdb):
    """Receptor의 중심, chain 중심, binding site 중심, Z범위 계산"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('receptor', receptor_pdb)
    
    # 전체 중심, chain별 중심, binding site 중심 계산
    all_coords = []
    chain_centers = {}
    binding_sites = {}
    
    for chain in structure.get_chains():
        chain_coords = []
        binding_coords = []
        
        for residue in chain.get_residues():
            # CA 원자 수집
            if 'CA' in residue:
                coord = residue['CA'].get_coord()
                chain_coords.append(coord)
                all_coords.append(coord)
                
                # Binding site 잔기인지 확인
                res_id = residue.get_id()[1]
                if res_id in BINDING_SITE_RESIDUES:
                    binding_coords.append(coord)
        
        if chain_coords:
            chain_centers[chain.id] = np.array(chain_coords).mean(axis=0)
        if binding_coords:
            binding_sites[chain.id] = np.array(binding_coords).mean(axis=0)
    
    all_coords = np.array(all_coords)
    total_center = all_coords.mean(axis=0)
    z_max = all_coords[:,2].max()  # Trimer의 최상단 Z 좌표
    z_min = all_coords[:,2].min()  # Trimer의 최하단 Z 좌표
    
    return total_center, chain_centers, binding_sites, z_max, z_min


def find_ligand_arms(mol, conf):
    """리간드의 3개 arm 끝점(L-glucose 위치) 찾기"""
    # 리간드 중심 계산
    coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    center = coords.mean(axis=0)
    
    # 중심에서 각 원자까지의 거리 계산
    distances = np.linalg.norm(coords - center, axis=1)
    
    # 산소 원자 중 가장 먼 3개 찾기 (L-glucose의 hydroxyl groups)
    o_atoms = []
    for i, atom in enumerate(mol.GetAtoms()):
        if atom.GetAtomicNum() == 8:  # Oxygen
            o_atoms.append((i, distances[i], coords[i]))
    
    # 거리순 정렬 후 가장 먼 원자들 선택
    o_atoms.sort(key=lambda x: x[1], reverse=True)
    
    # 3개의 arm 방향 벡터 찾기 (서로 다른 방향의 원자 선택)
    arm_points = []
    for idx, dist, coord in o_atoms:
        if len(arm_points) == 0:
            arm_points.append(coord)
        else:
            # 기존 arm들과 충분히 다른 방향인지 확인
            is_different = True
            for existing in arm_points:
                vec1 = (coord - center) / np.linalg.norm(coord - center)
                vec2 = (existing - center) / np.linalg.norm(existing - center)
                # 코사인 유사도가 0.5 미만이면 다른 방향
                if np.dot(vec1, vec2) > 0.5:
                    is_different = False
                    break
            if is_different:
                arm_points.append(coord)
        if len(arm_points) == 3:
            break
    
    return np.array(arm_points), center


def compute_rotation_matrix(ligand_arms, ligand_center, target_sites, total_center):
    """리간드 arm이 각 binding site를 향하도록 회전 행렬 계산"""
    # Binding site 방향 벡터 (중심에서 각 binding site로)
    site_dirs = []
    for chain_id in sorted(target_sites.keys()):
        direction = target_sites[chain_id] - total_center
        direction = direction / np.linalg.norm(direction)
        site_dirs.append(direction)
    chain_dirs = np.array(site_dirs)
    
    # 리간드 arm 방향 벡터
    arm_dirs = []
    for arm in ligand_arms:
        direction = arm - ligand_center
        direction = direction / np.linalg.norm(direction)
        arm_dirs.append(direction)
    arm_dirs = np.array(arm_dirs)
    
    # Kabsch 알고리즘으로 최적 회전 행렬 계산
    H = arm_dirs.T @ chain_dirs
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    
    # Reflection 보정
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    
    return R

def manual_docking(receptor_pdb, ligand_smiles, output_dir):
    print(f"📍 [Step 2 (Manual)] Placing Tripod manually...")
    
    # 1. 리간드 3D 구조 생성 (큰 분자용 설정)
    mol = Chem.MolFromSmiles(ligand_smiles)
    mol = Chem.AddHs(mol)
    
    # 큰 분자를 위한 EmbedMolecule 파라미터 설정
    params = AllChem.ETKDGv3()
    params.maxIterations = 5000
    params.randomSeed = 42
    params.useRandomCoords = True  # 랜덤 좌표로 시작 (큰 분자에 효과적)
    
    result = AllChem.EmbedMolecule(mol, params)
    if result == -1:
        print("   ⚠️ First embedding failed, trying with more iterations...")
        params.maxIterations = 10000
        params.useRandomCoords = True
        result = AllChem.EmbedMolecule(mol, params)
        if result == -1:
            raise ValueError("Failed to generate 3D coordinates for ligand")
    
    # 구조 최적화 (MMFF가 실패하면 UFF 사용)
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    except:
        print("   ⚠️ MMFF failed, using UFF force field...")
        AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
    
    # 2. Receptor 정보 계산
    total_center, chain_centers, binding_sites, z_max, z_min = get_receptor_info(receptor_pdb)
    print(f"   📍 Receptor center: {total_center}")
    print(f"   📍 Z range: {z_min:.1f} ~ {z_max:.1f} Å (막 관통 방향)")
    print(f"   📍 Binding sites found: {list(binding_sites.keys())}")
    
    # Binding site Z 좌표 확인 (채널 입구 = Z+ 방향)
    for chain_id, site in sorted(binding_sites.items()):
        print(f"      Chain {chain_id} binding site Z: {site[2]:.1f} Å")
    
    # 3. 리간드 arm 위치 찾기
    conf = mol.GetConformer()
    ligand_arms, ligand_center = find_ligand_arms(mol, conf)
    print(f"   📍 Found {len(ligand_arms)} ligand arms")
    
    # 4. 리간드를 먼저 중심으로 이동
    ligand_coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    ligand_center = ligand_coords.mean(axis=0)
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        new_pos = (pos.x - ligand_center[0], pos.y - ligand_center[1], pos.z - ligand_center[2])
        conf.SetAtomPosition(i, new_pos)
    
    # 5. Tripod를 GLUT1 trimer 위에 배치 (Z+ 방향)
    # Tripod core가 위에, L-glucose arm이 아래로 내려가도록
    # 리간드의 arm 방향을 아래(-Z)로 향하게 회전
    
    # 리간드 arm 방향 확인 및 회전
    ligand_coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    ligand_arms_new, _ = find_ligand_arms(mol, conf)
    
    # Arm들의 평균 Z 방향 확인
    arm_z_avg = ligand_arms_new[:,2].mean() if len(ligand_arms_new) > 0 else 0
    print(f"   📍 Ligand arms average Z: {arm_z_avg:.1f} (should be negative for downward)")
    
    # Arm이 위를 향하면 180도 회전 (X축 기준)
    if arm_z_avg > 0:
        print(f"   🔄 Flipping ligand (arms pointing up -> down)")
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            # X축 기준 180도 회전: (x, y, z) -> (x, -y, -z)
            conf.SetAtomPosition(i, (pos.x, -pos.y, -pos.z))
    
    # 6. Tripod를 trimer 위에 배치
    # 새 SMILES 구조: arm 길이 약 15Å (PEG6 + L-glucose)
    # Tripod core를 trimer 상단 바로 위에 배치
    binding_site_z = np.mean([site[2] for site in binding_sites.values()])
    tripod_height = z_max + 5  # Trimer 상단에서 5Å 위 (arm이 채널로 들어갈 수 있도록)
    
    print(f"   📍 Binding site Z: {binding_site_z:.1f} Å")
    print(f"   📍 Trimer Z_max: {z_max:.1f} Å")
    print(f"   📍 Target Tripod height: {tripod_height:.1f} Å")
    
    # XY 중심은 trimer 중심과 동일, Z는 계산된 높이
    target_pos = np.array([total_center[0], total_center[1], tripod_height])
    
    ligand_coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    current_center = ligand_coords.mean(axis=0)
    translation = target_pos - current_center
    
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        new_pos = (pos.x + translation[0], pos.y + translation[1], pos.z + translation[2])
        conf.SetAtomPosition(i, new_pos)
    
    # 최종 위치 확인
    final_coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    print(f"   📍 Tripod placed at Z={final_coords[:,2].mean():.1f} Å (above trimer Z_max={z_max:.1f})")
    print(f"   📍 Tripod Z range: {final_coords[:,2].min():.1f} ~ {final_coords[:,2].max():.1f} Å")
    
    # 6. 결과 폴더 생성 및 저장
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "rank1.sdf")
    
    w = Chem.SDWriter(output_file)
    w.write(mol)
    w.close()
    
    print(f"   ✅ Manual placement done: {output_file}")

if __name__ == "__main__":
    manual_docking(sys.argv[1], sys.argv[2], sys.argv[3])
