#!/usr/bin/env python3
"""
유연 폴리머를 위한 RDKit 3D 구조 생성

PEG24와 같은 초고유연 분자를 위한 임베딩 전략:
1. Heavy atom만으로 먼저 임베딩
2. useRandomCoords + EmbedMultipleConfs
3. 여러 컨포머 생성 후 최저 에너지 선택
4. 수소는 나중에 추가
"""

from pathlib import Path

def main():
    print("=" * 80)
    print("유연 폴리머 3D 구조 생성 (RDKit)")
    print("=" * 80)
    print()
    
    # 경로 설정
    data_dir = Path(__file__).parent.parent / "data"
    smiles_file = data_dir / "1arm_peg24_glc.smi"
    
    # SMILES 읽기
    print("Step 1: SMILES 읽기")
    print("-" * 80)
    with open(smiles_file, 'r') as f:
        smiles = f.read().strip()
    
    print(f"SMILES: {smiles[:80]}...")
    print(f"길이: {len(smiles)} 문자")
    print()
    
    # RDKit 임베딩
    print("Step 2: RDKit 유연 폴리머 임베딩")
    print("-" * 80)
    
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors
        
        # SMILES → Mol (수소 없이)
        m = Chem.MolFromSmiles(smiles)
        
        if m is None:
            print("❌ 유효하지 않은 SMILES")
            return
        
        print(f"✅ Mol 객체 생성 (heavy atoms: {m.GetNumAtoms()})")
        print()
        
        # 임베딩 파라미터 (유연 폴리머용)
        print("Step 3: 임베딩 파라미터 설정")
        print("-" * 80)
        
        params = AllChem.ETKDGv3()
        params.randomSeed = 0xF00D
        params.useRandomCoords = True
        params.pruneRmsThresh = 0.5
        params.numThreads = 0  # 모든 코어 사용
        
        print("파라미터:")
        print(f"  - useRandomCoords: {params.useRandomCoords}")
        print(f"  - pruneRmsThresh: {params.pruneRmsThresh}")
        print()
        
        # 여러 컨포머 생성
        print("Step 4: 다중 컨포머 임베딩 (시간이 걸립니다)")
        print("-" * 80)
        print("50개 컨포머 시도 중...")
        
        cids = AllChem.EmbedMultipleConfs(m, numConfs=50, params=params)
        
        if not cids:
            print("❌ 모든 컨포머 임베딩 실패")
            print("   대안: OpenBabel 사용 또는 파라미터 조정")
            return
        
        print(f"✅ {len(cids)}개 컨포머 생성 성공")
        print()
        
        # UFF 최적화
        print("Step 5: UFF 최적화")
        print("-" * 80)
        
        success_count = 0
        for i, cid in enumerate(cids):
            try:
                result = AllChem.UFFOptimizeMolecule(m, confId=cid, maxIters=5000)
                if result == 0:
                    success_count += 1
                if (i + 1) % 10 == 0:
                    print(f"  진행: {i+1}/{len(cids)} 완료")
            except Exception as e:
                pass
        
        print(f"✅ {success_count}/{len(cids)} 컨포머 최적화 성공")
        print()
        
        # 최저 에너지 컨포머 선택
        print("Step 6: 최저 에너지 컨포머 선택")
        print("-" * 80)
        
        energies = []
        for cid in cids:
            try:
                ff = AllChem.UFFGetMoleculeForceField(m, confId=cid)
                e = ff.CalcEnergy()
                energies.append((e, cid))
            except Exception:
                energies.append((1e18, cid))
        
        energies.sort()
        best_energy, best_cid = energies[0]
        
        print(f"최저 에너지: {best_energy:.2f} kcal/mol (conformer {best_cid})")
        print(f"최고 에너지: {energies[-1][0]:.2f} kcal/mol")
        print()
        
        # 최적 컨포머 추출
        print("Step 7: 최적 컨포머 추출 및 수소 추가")
        print("-" * 80)
        
        best = Chem.Mol(m)
        best.RemoveAllConformers()
        best.AddConformer(m.GetConformer(best_cid), assignId=True)
        
        # 이제 수소 추가
        best_h = Chem.AddHs(best, addCoords=True)
        
        print(f"✅ 수소 추가 완료 (총 원자: {best_h.GetNumAtoms()})")
        print()
        
        # 파일 저장
        print("Step 8: 파일 저장")
        print("-" * 80)
        
        # MOL2 (CGenFF 입력용)
        mol2_file = data_dir / "1arm_peg24_glc.mol2"
        Chem.MolToMolFile(best_h, str(mol2_file))
        print(f"✅ MOL2: {mol2_file}")
        
        # PDB (시각화용)
        pdb_file = data_dir / "1arm_peg24_glc.pdb"
        Chem.MolToPDBFile(best_h, str(pdb_file))
        print(f"✅ PDB: {pdb_file}")
        
        # SDF (백업)
        sdf_file = data_dir / "1arm_peg24_glc.sdf"
        writer = Chem.SDWriter(str(sdf_file))
        writer.write(best_h)
        writer.close()
        print(f"✅ SDF: {sdf_file}")
        print()
        
        # 분자 정보
        print("Step 9: 최종 분자 정보")
        print("-" * 80)
        print(f"  - 분자식: {Chem.rdMolDescriptors.CalcMolFormula(best_h)}")
        print(f"  - 분자량: {Descriptors.MolWt(best_h):.2f} g/mol")
        print(f"  - 원자 수: {best_h.GetNumAtoms()}")
        print(f"  - 회전 가능 결합: {Descriptors.NumRotatableBonds(best_h)}")
        print()
        
    except ImportError:
        print("❌ RDKit이 설치되지 않았습니다.")
        return
    
    # 완료
    print("=" * 80)
    print("✅ 3D 구조 생성 완료!")
    print("=" * 80)
    print()
    print("생성된 파일:")
    print(f"  - {mol2_file} (CGenFF 입력용)")
    print(f"  - {pdb_file} (시각화용)")
    print(f"  - {sdf_file} (백업)")
    print()
    print("다음 단계:")
    print("  1. 구조 확인:")
    print("     pymol data/1arm_peg24_glc.pdb")
    print()
    print("  2. CGenFF 파라미터 생성:")
    print("     https://cgenff.umaryland.edu/")
    print("     - MOL2 파일 업로드")
    print("     - .str 파일 다운로드")
    print()
    print("  3. Penalty score 확인:")
    print("     grep 'PENALTY' data/*.str | sort -k2 -n -r | head -20")
    print()


if __name__ == "__main__":
    main()
