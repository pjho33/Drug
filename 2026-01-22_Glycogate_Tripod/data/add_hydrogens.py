#!/usr/bin/env python3
"""
SMILES에 수소 추가
RDKit 사용
"""

from rdkit import Chem
from pathlib import Path

def add_hydrogens_to_smiles(input_smi, output_smi):
    """SMILES 파일을 읽어서 수소를 추가하고 저장"""
    
    # SMILES 읽기
    with open(input_smi, 'r') as f:
        smiles = f.read().strip()
    
    print(f"원본 SMILES: {smiles[:100]}...")
    
    # RDKit 분자 객체 생성
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        print("ERROR: SMILES 파싱 실패!")
        return False
    
    print(f"원자 수 (수소 제외): {mol.GetNumAtoms()}")
    
    # 수소 추가
    mol_with_h = Chem.AddHs(mol)
    print(f"원자 수 (수소 포함): {mol_with_h.GetNumAtoms()}")
    
    # 수소 포함 SMILES 생성
    smiles_with_h = Chem.MolToSmiles(mol_with_h)
    
    # 저장
    with open(output_smi, 'w') as f:
        f.write(smiles_with_h)
    
    print(f"\n수소 추가 완료!")
    print(f"출력 파일: {output_smi}")
    print(f"수소 포함 SMILES: {smiles_with_h[:100]}...")
    
    return True

if __name__ == "__main__":
    # 파일 경로
    data_dir = Path(__file__).parent
    input_file = data_dir / "3arm_peg24_glc.smi"
    output_file = data_dir / "3arm_peg24_glc_with_h.smi"
    
    print("=" * 80)
    print("3-Arm Tripod SMILES에 수소 추가")
    print("=" * 80)
    print()
    
    if not input_file.exists():
        print(f"ERROR: 입력 파일을 찾을 수 없습니다: {input_file}")
        exit(1)
    
    success = add_hydrogens_to_smiles(input_file, output_file)
    
    if success:
        print("\n✅ 성공!")
        print(f"\nCGenFF 업로드 파일: {output_file.name}")
    else:
        print("\n❌ 실패!")
        exit(1)
