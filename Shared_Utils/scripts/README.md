# GLUT1-Targeting Drug Design Project - Scripts Organization

이 디렉토리는 프로젝트 진행 과정에서 생성된 모든 스크립트를 **날짜와 목적별로 정리**한 것입니다.

## 📋 프로젝트 진행 순서

### Phase 1: 리간드 설계 및 도킹 (2025-12)
1. **01_Phase1_Ligand_Design** - Short-Tripod 리간드 SMILES 생성
2. **02_Phase1_Docking** - GLUT1 단백질 준비 및 AutoDock Vina 도킹

### Phase 2: MD 시뮬레이션 및 분석 (2025-12 ~ 2026-01)
3. **03_Phase2_MD_Simulation** - OpenMM을 이용한 MD 시뮬레이션
4. **04_Phase2_Analysis** - RMSD, 접촉 분석, MM/GBSA 계산

### Phase 3: CHARMM-GUI 복합체 준비 (2026-01)
5. **05_Complex_Preparation** - GLUT1-Ligand 복합체 PDB 생성
6. **06_PDB_Format_Fixing** - CHARMM-GUI 호환 포맷 수정
7. **07_PSF_Generation** - CHARMM topology 파일 생성

### Phase 4: Control MD 및 분석 (2026-01-08)
8. **08_Control_MD_CHARMM** - CHARMM-GUI OpenMM으로 100ns MD
9. **09_Control_MD_Analysis** - RMSD, RMSF, 수소결합 분석

### Phase 5: MM/PBSA 결합 에너지 (2026-01-10~11)
10. **10_MMPBSA_Attempts** - AmberTools MMPBSA.py를 이용한 결합 에너지 계산

## 📊 주요 결과

### Phase 2 결과
- 리간드 안정성 확인
- RMSD < 2Å 유지
- MM/GBSA 결합 에너지 계산

### Control MD 결과 (100ns)
- Protein Backbone RMSD: 1.36 ± 0.52 Å
- Ligand RMSD: 1.18 ± 0.42 Å
- Mean RMSF: 1.17 ± 0.32 Å
- Ligand Burial: 99.3%
- 추정 결합 에너지: -576.50 kcal/mol

## 🔧 기술 스택

- **분자 모델링:** RDKit, OpenBabel
- **도킹:** AutoDock Vina
- **MD 시뮬레이션:** OpenMM (CUDA), CHARMM-GUI
- **Force Field:** CHARMM36, GAFF2
- **분석:** MDAnalysis, PyMOL
- **결합 에너지:** AmberTools (MMPBSA.py)

## 📝 참고사항

- 원본 스크립트는 `/home/pjho3/projects/Drug/scripts/`에 그대로 유지됩니다
- 이 정리된 폴더는 **참조 및 문서화 목적**입니다
- 각 폴더의 README.md에 상세 설명이 있습니다

---

**Last Updated:** 2026-01-11
