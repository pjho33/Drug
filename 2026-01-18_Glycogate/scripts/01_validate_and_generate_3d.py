#!/usr/bin/env python3
"""
1-Arm SMILES 검증 및 3D 구조 생성

TRIS-PEG24-β‑L‑Glc 분자의 SMILES를 검증하고
3D 구조를 생성합니다.
"""

import os
from pathlib import Path

def main():
    print("=" * 80)
    print("1-Arm SMILES 검증 및 3D 구조 생성")
    print("=" * 80)
    print()
    
    # 경로 설정
    data_dir = Path(__file__).parent.parent / "data"
    smiles_file = data_dir / "1arm_smiles.txt"
    
    # SMILES 읽기
    print("Step 1: SMILES 읽기")
    print("-" * 80)
    with open(smiles_file, 'r') as f:
        smiles = f.read().strip()
    
    print(f"SMILES: {smiles}")
    print(f"길이: {len(smiles)} 문자")
    print()
    
    # RDKit으로 검증
    print("Step 2: RDKit 검증")
    print("-" * 80)
    
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors
        
        # SMILES → Mol 객체
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            print("❌ 오류: 유효하지 않은 SMILES")
            return
        
        print("✅ 유효한 SMILES")
        print()
        
        # 분자 정보
        print("분자 정보:")
        print(f"  - 분자식: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
        print(f"  - 분자량: {Descriptors.MolWt(mol):.2f} g/mol")
        print(f"  - 원자 수: {mol.GetNumAtoms()}")
        print(f"  - 결합 수: {mol.GetNumBonds()}")
        print(f"  - 회전 가능 결합: {Descriptors.NumRotatableBonds(mol)}")
        print(f"  - H-bond donor: {Descriptors.NumHDonors(mol)}")
        print(f"  - H-bond acceptor: {Descriptors.NumHAcceptors(mol)}")
        print()
        
        # 원자 종류별 개수
        atom_counts = {}
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
        
        print("원자 구성:")
        for symbol, count in sorted(atom_counts.items()):
            print(f"  - {symbol}: {count}")
        print()
        
    except ImportError:
        print("⚠️  RDKit이 설치되지 않았습니다.")
        print("   conda install -c conda-forge rdkit")
        print()
        return
    
    # 3D 구조 생성
    print("Step 3: 3D 구조 생성")
    print("-" * 80)
    
    # 수소 추가
    mol_h = Chem.AddHs(mol)
    print(f"✅ 수소 추가 완료 (총 원자 수: {mol_h.GetNumAtoms()})")
    
    # 3D 좌표 생성 (ETKDG)
    print("3D 좌표 생성 중 (ETKDG 알고리즘)...")
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.numThreads = 0  # 모든 코어 사용
    
    result = AllChem.EmbedMolecule(mol_h, params)
    
    if result == -1:
        print("❌ 3D 좌표 생성 실패")
        print("   분자가 너무 크거나 복잡할 수 있습니다.")
        print("   대안: 단계적 구조 생성 또는 외부 도구 사용")
        return
    
    print("✅ 3D 좌표 생성 완료")
    print()
    
    # 구조 최적화 (UFF)
    print("Step 4: 구조 최적화 (UFF)")
    print("-" * 80)
    print("최적화 중 (시간이 걸릴 수 있습니다)...")
    
    try:
        result = AllChem.UFFOptimizeMolecule(mol_h, maxIters=2000)
        print(f"✅ 최적화 완료 (수렴 코드: {result})")
        print("   0 = 수렴, 1 = 최대 반복 도달")
        print()
    except Exception as e:
        print(f"⚠️  최적화 중 오류: {e}")
        print("   최적화 없이 진행합니다.")
        print()
    
    # 파일 저장
    print("Step 5: 파일 저장")
    print("-" * 80)
    
    # MOL2 파일 (CGenFF 입력용)
    mol2_file = data_dir / "1arm_peg24_glc.mol2"
    Chem.MolToMolFile(mol_h, str(mol2_file))
    print(f"✅ MOL2 파일: {mol2_file}")
    
    # PDB 파일 (시각화용)
    pdb_file = data_dir / "1arm_peg24_glc.pdb"
    Chem.MolToPDBFile(mol_h, str(pdb_file))
    print(f"✅ PDB 파일: {pdb_file}")
    
    # SDF 파일 (백업)
    sdf_file = data_dir / "1arm_peg24_glc.sdf"
    writer = Chem.SDWriter(str(sdf_file))
    writer.write(mol_h)
    writer.close()
    print(f"✅ SDF 파일: {sdf_file}")
    print()
    
    # 완료
    print("=" * 80)
    print("✅ 완료!")
    print("=" * 80)
    print()
    print("다음 단계:")
    print("  1. 구조 확인: pymol 1arm_peg24_glc.pdb")
    print("  2. CGenFF 파라미터 생성:")
    print("     - CGenFF Server: https://cgenff.umaryland.edu/")
    print("     - 또는 로컬 CGenFF 프로그램 사용")
    print("  3. Penalty score 확인")
    print()


if __name__ == "__main__":
    main()
