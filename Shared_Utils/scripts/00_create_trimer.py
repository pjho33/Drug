# scripts/00_create_trimer.py
import sys
import math
import copy
import numpy as np
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model

def rotate_chain_x90(chain):
    """X축 기준 -90도 회전 (Y→Z, Z→-Y) - 막 관통 방향을 Z축으로"""
    for atom in chain.get_atoms():
        x, y, z = atom.coord
        # X축 기준 -90도 회전: (x, y, z) → (x, z, -y)
        atom.coord = np.array([x, z, -y])

def center_chain(chain):
    """체인을 원점 중심으로 이동"""
    coords = np.array([atom.coord for atom in chain.get_atoms()])
    center = coords.mean(axis=0)
    for atom in chain.get_atoms():
        atom.coord = atom.coord - center

def create_trimer(input_pdb, output_pdb, distance_angstrom):
    print(f"🏗️ [Modeling] Creating Trimer from {input_pdb}...")
    print(f"   📐 Distance between units: {distance_angstrom} Å")

    # 1. 원본 파일 로딩
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('monomer', input_pdb)
    model = structure[0]

    # 2. 새로운 구조체 생성
    new_structure = Structure("Trimer")
    new_model = Model(0)
    new_structure.add(new_model)

    # 3. 배치 좌표 계산 (XY 평면에 정삼각형, 위에서 보면 삼각형)
    # 막 관통 방향 = Z축 (glucose가 위에서 아래로 통과)
    d = distance_angstrom
    # 정삼각형 중심이 원점에 오도록 배치
    positions = [
        (-d/2.0, -d * math.sqrt(3) / 6.0, 0.0),  # 왼쪽 아래
        (d/2.0, -d * math.sqrt(3) / 6.0, 0.0),   # 오른쪽 아래
        (0.0, d * math.sqrt(3) / 3.0, 0.0)       # 위쪽
    ]
    
    chain_ids = ['A', 'B', 'C']

    # 원본 체인 가져오기
    original_chain = list(model.get_chains())[0]

    # 4. 복제, 회전, 이동
    for i, (x_shift, y_shift, z_shift) in enumerate(positions):
        print(f"   👉 Generating Unit {i+1} (Chain {chain_ids[i]})...")
        
        # 체인 복사
        new_chain = copy.deepcopy(original_chain)
        new_chain.id = chain_ids[i]
        
        # 먼저 중심을 원점으로
        center_chain(new_chain)
        
        # X축 기준 -90도 회전 (Y축이 막 관통 방향 → Z축이 막 관통 방향)
        rotate_chain_x90(new_chain)
        
        # 정삼각형 위치로 이동
        for atom in new_chain.get_atoms():
            atom.coord[0] += x_shift
            atom.coord[1] += y_shift
            atom.coord[2] += z_shift
            
        new_model.add(new_chain)
    
    # 회전 후 좌표 범위 확인
    all_coords = np.array([atom.coord for atom in new_model.get_atoms()])
    print(f"   📊 Trimer 좌표 범위:")
    print(f"      X: {all_coords[:,0].min():.1f} ~ {all_coords[:,0].max():.1f} Å")
    print(f"      Y: {all_coords[:,1].min():.1f} ~ {all_coords[:,1].max():.1f} Å")
    print(f"      Z: {all_coords[:,2].min():.1f} ~ {all_coords[:,2].max():.1f} Å (막 관통 방향)")

    # 5. 저장
    io = PDBIO()
    io.set_structure(new_structure)
    io.save(output_pdb)
    print("-" * 50)
    print(f"   ✅ Trimer PDB saved to: {output_pdb}")
    print("-" * 50)

if __name__ == "__main__":
    # 사용법: python 00_create_trimer.py <입력파일> <출력파일> <거리>
    if len(sys.argv) != 4:
        print("Usage: python 00_create_trimer.py <input.pdb> <output.pdb> <distance>")
        sys.exit(1)

    create_trimer(sys.argv[1], sys.argv[2], float(sys.argv[3]))