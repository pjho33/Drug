#!/usr/bin/env python3
"""
3-arm에서 1-arm SMILES 추출

원본 Tris-(PEG2-Glc)3에서 1개 팔만 남기고 나머지는 OH로 캡핑
"""

from pathlib import Path

def main():
    print("=" * 80)
    print("3-arm → 1-arm SMILES 변환")
    print("=" * 80)
    print()
    
    # 원본 3-arm SMILES
    original_3arm = "NC(COCCOCCOO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O)(COCCOCCOO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O)COCCOCCOO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O"
    
    print("Step 1: 원본 3-arm 구조 분석")
    print("-" * 80)
    print("구조: NC(arm1)(arm2)arm3")
    print()
    print("각 팔:")
    arm = "COCCOCCOO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O"
    print(f"  {arm}")
    print()
    
    # 1-arm SMILES 생성
    print("Step 2: 1-arm SMILES 생성")
    print("-" * 80)
    
    # TRIS 중심: NC(CO)(CO)C...
    # 1개 팔만 남기고 나머지 2개는 OH로 캡핑
    smiles_1arm = f"NC(CO)(CO){arm}"
    
    print("1-arm SMILES:")
    print(smiles_1arm)
    print()
    
    # RDKit으로 검증
    print("Step 3: RDKit 검증")
    print("-" * 80)
    
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw, Descriptors
        
        mol = Chem.MolFromSmiles(smiles_1arm)
        
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
        
        # 구조 이미지 저장
        data_dir = Path(__file__).parent.parent / "data"
        img_file = data_dir / "1arm_peg2_glc_structure.png"
        
        img = Draw.MolToImage(mol, size=(800, 600))
        img.save(str(img_file))
        print(f"✅ 구조 이미지: {img_file}")
        print()
        
        # SMILES 파일 저장
        smiles_file = data_dir / "1arm_peg2_glc.smi"
        with open(smiles_file, 'w') as f:
            f.write(smiles_1arm)
        print(f"✅ SMILES 저장: {smiles_file}")
        print()
        
    except ImportError:
        print("❌ RDKit이 설치되지 않았습니다.")
        return
    
    # 완료
    print("=" * 80)
    print("✅ 1-arm SMILES 생성 완료!")
    print("=" * 80)
    print()
    print("현재 구조: TRIS(OH)2 - PEG2 - O - Glucose")
    print()
    print("다음 단계:")
    print("  1. 구조 확인: 이미지 파일 열기")
    print("  2. PEG2 → PEG24 확장")
    print("  3. Triazole-Amide 링커 추가 (필요시)")
    print()


if __name__ == "__main__":
    main()
