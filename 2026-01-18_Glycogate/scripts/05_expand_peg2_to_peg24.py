#!/usr/bin/env python3
"""
PEG2를 PEG24로 확장

1-arm TRIS(OH)2-PEG2-Glc를 TRIS(OH)2-PEG24-Glc로 확장합니다.
"""

from pathlib import Path

def main():
    print("=" * 80)
    print("PEG2 → PEG24 확장")
    print("=" * 80)
    print()
    
    # 기존 1-arm PEG2 SMILES 읽기
    data_dir = Path(__file__).parent.parent / "data"
    peg2_file = data_dir / "1arm_peg2_glc.smi"
    
    with open(peg2_file, 'r') as f:
        smiles_peg2 = f.read().strip()
    
    print("Step 1: 기존 PEG2 SMILES")
    print("-" * 80)
    print(smiles_peg2)
    print()
    
    # PEG2 부분 분석
    print("Step 2: PEG 부분 분석")
    print("-" * 80)
    print("현재 PEG2: COCCOC")
    print("  - C-O-C-C-O-C")
    print("  - 2개 EO 유닛")
    print()
    
    # PEG24 생성
    print("Step 3: PEG24 생성")
    print("-" * 80)
    
    # PEG 유닛: OCCOC (한 개 EO = -O-CH2-CH2-O-C)
    # 24개 반복
    peg24_chain = "OCCOC" * 24
    
    print(f"PEG24 체인 길이: {len(peg24_chain)} 문자")
    print(f"PEG24: {peg24_chain[:50]}... (처음 50자)")
    print()
    
    # 1-arm PEG24 SMILES 조합
    print("Step 4: 1-arm PEG24 SMILES 조합")
    print("-" * 80)
    
    # TRIS 중심
    tris_core = "NC(CO)(CO)C"
    
    # Glucose (원본에서 추출)
    glucose = "O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O"
    
    # 조합: TRIS - PEG24 - O - Glucose
    smiles_peg24 = f"{tris_core}{peg24_chain}O{glucose}"
    
    print("생성된 SMILES:")
    print(smiles_peg24)
    print()
    print(f"총 길이: {len(smiles_peg24)} 문자")
    print()
    
    # RDKit으로 검증
    print("Step 5: RDKit 검증")
    print("-" * 80)
    
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw, Descriptors
        
        mol = Chem.MolFromSmiles(smiles_peg24)
        
        if mol is None:
            print("❌ 유효하지 않은 SMILES")
            return
        
        print("✅ 유효한 SMILES")
        print()
        
        # 분자 정보
        print("분자 정보:")
        print(f"  - 분자식: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
        print(f"  - 분자량: {Descriptors.MolWt(mol):.2f} g/mol")
        print(f"  - 원자 수: {mol.GetNumAtoms()}")
        print(f"  - 회전 가능 결합: {Descriptors.NumRotatableBonds(mol)}")
        print(f"  - H-bond donor: {Descriptors.NumHDonors(mol)}")
        print(f"  - H-bond acceptor: {Descriptors.NumHAcceptors(mol)}")
        print()
        
        # 원자 구성
        atom_counts = {}
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
        
        print("원자 구성:")
        for symbol, count in sorted(atom_counts.items()):
            print(f"  - {symbol}: {count}")
        print()
        
        # 예상 원자 수 확인
        print("예상 vs 실제:")
        print(f"  - C: 예상 ~63 (TRIS 3 + PEG 48 + Glc 6), 실제 {atom_counts.get('C', 0)}")
        print(f"  - O: 예상 ~33 (TRIS 2 + PEG 25 + Glc 6), 실제 {atom_counts.get('O', 0)}")
        print(f"  - N: 예상 1, 실제 {atom_counts.get('N', 0)}")
        print()
        
        # SMILES 파일 저장
        smiles_file = data_dir / "1arm_peg24_glc.smi"
        with open(smiles_file, 'w') as f:
            f.write(smiles_peg24)
        print(f"✅ SMILES 저장: {smiles_file}")
        print()
        
        # 구조 이미지는 너무 커서 생략 (261 원자)
        print("⚠️  구조 이미지 생성 생략 (분자가 너무 큼)")
        print("   3D 구조 생성 후 PyMOL/VMD로 확인하세요.")
        print()
        
    except ImportError:
        print("❌ RDKit이 설치되지 않았습니다.")
        return
    
    # 완료
    print("=" * 80)
    print("✅ PEG24 확장 완료!")
    print("=" * 80)
    print()
    print("현재 구조: TRIS(OH)2 - PEG24 - O - Glucose")
    print()
    print("다음 단계:")
    print("  1. 3D 구조 생성 (OpenBabel 사용)")
    print("  2. CGenFF 파라미터 생성")
    print("  3. MD 시뮬레이션 준비")
    print()
    print("참고:")
    print("  - 이 구조는 Triazole-Amide 링커가 없는 단순 PEG-O-Glc입니다.")
    print("  - 실제 합성 구조와 다를 수 있으니 확인이 필요합니다.")
    print()


if __name__ == "__main__":
    main()
