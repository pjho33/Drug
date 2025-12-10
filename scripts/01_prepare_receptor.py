# scripts/01_prepare_receptor.py
import sys
import os
from Bio import PDB
from pdbfixer import PDBFixer
from openmm.app import PDBFile

class ChainSelect(PDB.Select):
    def accept_chain(self, chain):
        # 첫 번째 체인(Chain A)만 선택하여 남깁니다.
        # 복잡한 Dimer/Trimer 구조에서 에러를 방지하는 핵심 로직입니다.
        model = chain.get_parent()
        first_chain_id = list(model.get_chains())[0].id
        return chain.id == first_chain_id

def prepare_receptor(input_pdb, output_pdb):
    print(f"🛠️ [Step 1] Preparing Receptor (High-End Mode): {input_pdb}")
    
    # ----------------------------------------------------
    # 1. Biopython으로 Chain A 안전 추출 (Extraction)
    # ----------------------------------------------------
    print("   ✂️ Extracting Chain A using Biopython...")
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)
    
    io = PDB.PDBIO()
    io.set_structure(structure)
    
    temp_pdb = "temp_chain_a.pdb"
    io.save(temp_pdb, ChainSelect())
    
    # ----------------------------------------------------
    # 2. PDBFixer로 수리 (Fixing)
    # ----------------------------------------------------
    print("   🔧 Loading extracted chain into PDBFixer...")
    fixer = PDBFixer(filename=temp_pdb)
    
    # 청소: 물, 이온, 잡동사니 제거
    print("   🧹 Removing artifacts (Water, Ions, Heterogens)...")
    fixer.removeHeterogens(keepWater=False)
    
    # 수리: 빠진 잔기 및 원자 복구
    print("   🔧 Fixing missing residues and atoms...")
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    
    # 수소 추가: pH 7.0 기준 (표준 생체 조건)
    print("   💧 Adding Hydrogens (pH 7.0)...")
    fixer.addMissingHydrogens(7.0)
    
    # ----------------------------------------------------
    # 3. 최종 저장
    # ----------------------------------------------------
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    
    print(f"   ✅ Saved cleaned receptor to: {output_pdb}")
    
    # 임시 파일 정리
    if os.path.exists(temp_pdb):
        os.remove(temp_pdb)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python 01_prepare_receptor.py <input_pdb> <output_pdb>")
        sys.exit(1)
    
    prepare_receptor(sys.argv[1], sys.argv[2])