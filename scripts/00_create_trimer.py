# scripts/00_create_trimer.py
import sys
import math
import copy
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model

def create_trimer(input_pdb, output_pdb, distance_angstrom):
    print(f"🏗️ [Modeling] Creating Trimer from {input_pdb}...")
    print(f"   📐 Distance between units: {distance_angstrom} Å")

    # 1. 원본 파일 로딩
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('monomer', input_pdb)
    model = structure[0] # 첫 번째 모델 사용

    # 2. 새로운 구조체 생성 (여기에 3개를 담을 예정)
    new_structure = Structure("Trimer")
    new_model = Model(0)
    new_structure.add(new_model)

    # 3. 배치 좌표 계산 (정삼각형)
    # Unit 1: 원점 (0, 0, 0)
    # Unit 2: X축으로 distance만큼 이동
    # Unit 3: 위쪽으로 이동 (정삼각형 꼭지점)
    d = distance_angstrom
    positions = [
        (0.0, 0.0, 0.0),
        (d, 0.0, 0.0),
        (d / 2.0, d * math.sqrt(3) / 2.0, 0.0)
    ]
    
    chain_ids = ['A', 'B', 'C'] # 체인 이름 변경 (A, B, C)

    # 원본 체인 가져오기 (첫 번째 체인만 사용)
    original_chain = list(model.get_chains())[0]

    # 4. 복제 및 이동
    for i, (x_shift, y_shift, z_shift) in enumerate(positions):
        print(f"   👉 Generating Unit {i+1} (Chain {chain_ids[i]})...")
        
        # 체인 복사 (Deep Copy: 원본 훼손 방지)
        new_chain = copy.deepcopy(original_chain)
        new_chain.id = chain_ids[i] # 체인 ID 변경
        
        # 원자 좌표 이동 (Translation)
        for atom in new_chain.get_atoms():
            atom.coord[0] += x_shift
            atom.coord[1] += y_shift
            atom.coord[2] += z_shift
            
        # 새 모델에 추가
        new_model.add(new_chain)

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