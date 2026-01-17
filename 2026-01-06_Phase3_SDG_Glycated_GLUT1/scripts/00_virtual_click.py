# scripts/00_virtual_click.py
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

def perform_virtual_click(core_smiles, arm_smiles, output_file):
    print(f"⚗️ [Virtual Synthesis] CuAAC Click Reaction Starting...")
    print(f"   Core (Azide): {core_smiles}")
    print(f"   Arm (Alkyne): {arm_smiles}")

    # 1. 분자 객체 생성
    core = Chem.MolFromSmiles(core_smiles)
    arm = Chem.MolFromSmiles(arm_smiles)

    # 2. 반응 정의 (CuAAC: Azide + Terminal Alkyne -> 1,4-Triazole)
    # SMARTS 패턴: Azide([N]=[N+]=[N-]) + Alkyne(C#C) -> Triazole
    # 이 패턴은 1,4-disubstituted 1,2,3-triazole을 형성합니다.
    rxn_smarts = '[*:1][N]=[N+]=[N-].[*:2][C]#[C]>>[*:1]n1cc([*:2])nn1'
    rxn = AllChem.ReactionFromSmarts(rxn_smarts)

    # 3. 반응 수행 (반복적으로 모든 자리에 붙임)
    # 코어에 Azide가 여러 개면 순차적으로 붙입니다.
    product = core
    
    # 코어에 있는 Azide 개수 확인
    azide_pattern = Chem.MolFromSmarts('[N]=[N+]=[N-]')
    num_sites = len(core.GetSubstructMatches(azide_pattern))
    print(f"   👉 Found {num_sites} conjugation sites (Azides) on Core.")

    current_mol = core
    for i in range(num_sites):
        # 반응 실행
        results = rxn.RunReactants((current_mol, arm))
        if not results:
            print("   ⚠️ Reaction failed at step", i+1)
            break
        
        # 첫 번째 생성물을 선택 (보통 이성질체 중 하나)
        current_mol = results[0][0]
        # 화학적 오류 수정 (Sanitize)
        Chem.SanitizeMol(current_mol)
        print(f"      🔗 Attached Arm {i+1}/{num_sites}...")

    # 4. 최종 결과 정리
    final_smiles = Chem.MolToSmiles(current_mol)
    print("-" * 50)
    print(f"   ✅ Final Tripod SMILES Generated!")
    print(f"   🧬 SMILES: {final_smiles}")
    print("-" * 50)

    # 파일로 저장 (SMILES 파일)
    with open(output_file, 'w') as f:
        f.write(final_smiles)
    print(f"   💾 Saved to: {output_file}")

if __name__ == "__main__":
    # 사용법: python 00_virtual_click.py "코어SMILES" "약물SMILES" "결과파일"
    if len(sys.argv) != 4:
        print("Usage: python 00_virtual_click.py <Core_Azide> <Arm_Alkyne> <Output.smi>")
        sys.exit(1)

    perform_virtual_click(sys.argv[1], sys.argv[2], sys.argv[3])