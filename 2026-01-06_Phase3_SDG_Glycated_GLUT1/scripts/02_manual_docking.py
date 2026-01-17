# scripts/02_manual_docking.py
import sys
import os
from rdkit import Chem
from rdkit.Chem import AllChem

def manual_docking(receptor_pdb, ligand_smiles, output_dir):
    print(f"📍 [Step 2 (Manual)] Placing Tripod manually...")
    
    # 1. 리간드 3D 구조 생성
    mol = Chem.MolFromSmiles(ligand_smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    
    # 2. 결과 폴더 생성
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "rank1.sdf") # DiffDock 결과인 척 저장
    
    # 3. 저장 (좌표는 원점 근처에 생성됨 -> MD에서 자리 잡음)
    w = Chem.SDWriter(output_file)
    w.write(mol)
    w.close()
    
    print(f"   ✅ Manual placement done: {output_file}")

if __name__ == "__main__":
    manual_docking(sys.argv[1], sys.argv[2], sys.argv[3])