# AnchorR 분석에 필요한 자료

## 필수 자료 (2개)

### 1. Topology 파일 (원자/결합/잔기 정보)

**파일 형식:**
- `.psf` (CHARMM-GUI/OpenMM) ← **가장 일반적**
- `.pdb` (단, 결합 정보 포함 필요)
- `.prmtop` (AMBER)

**위치 예시:**
```
~/projects/Drug/2026-01-22_Glycogate_Tripod/data/solution_builder/openmm/step3_input.psf
```

**확인 방법:**
```bash
ls -lh ~/projects/Drug/2026-01-22_Glycogate_Tripod/data/solution_builder/openmm/*.psf
```

---

### 2. Trajectory 파일 (시간에 따른 좌표)

**파일 형식:**
- `.dcd` (CHARMM/OpenMM) ← **가장 일반적**
- `.xtc` (GROMACS)
- `.nc` (AMBER)

**위치 예시:**
```
~/projects/Drug/2026-01-22_Glycogate_Tripod/results/md_3arm_tripod/step5_production.dcd
```

**확인 방법:**
```bash
ls -lh ~/projects/Drug/2026-01-22_Glycogate_Tripod/results/md_3arm_tripod/*.dcd
```

---

## 선택 자료 (있으면 좋음)

### 3. Reference 좌표 파일

**파일 형식:**
- `.pdb` (처음 프레임)
- `.inpcrd` / `.rst7` (AMBER)

**용도:**
- Topology가 `.prmtop`일 때 필요할 수 있음
- CHARMM/OpenMM에서는 보통 불필요

---

## 현재 상황별 준비 방법

### Case 1: Tripod MD 시뮬레이션이 **완료된 경우**

**필요한 파일:**
```bash
# Topology
~/projects/Drug/2026-01-22_Glycogate_Tripod/data/solution_builder/openmm/step3_input.psf

# Trajectory
~/projects/Drug/2026-01-22_Glycogate_Tripod/results/md_3arm_tripod/step5_production.dcd
```

**확인 명령:**
```bash
# PSF 파일 확인
ls -lh ~/projects/Drug/2026-01-22_Glycogate_Tripod/data/solution_builder/openmm/step3_input.psf

# DCD 파일 확인
ls -lh ~/projects/Drug/2026-01-22_Glycogate_Tripod/results/md_3arm_tripod/*.dcd

# 파일 크기 확인 (DCD는 보통 수 GB)
du -h ~/projects/Drug/2026-01-22_Glycogate_Tripod/results/md_3arm_tripod/*.dcd
```

---

### Case 2: Tripod MD 시뮬레이션이 **진행 중인 경우**

**옵션 1: 일부 trajectory로 테스트**
- 처음 10,000 프레임만 사용
- 분석 스크립트 테스트 및 파라미터 조정

**옵션 2: 시뮬레이션 완료 대기**
- 20ns 전체 완료 후 분석
- 더 정확한 통계

---

### Case 3: Tripod MD 시뮬레이션이 **아직 시작 안 된 경우**

**준비 순서:**
1. Tripod SMILES 생성 ✅ (완료)
2. CGenFF 파라미터 생성
3. Ligand Reader 실행
4. Solution Builder 실행
5. 20ns MD 시뮬레이션 실행
6. **분석 시작** ← 여기서 AnchorR 스크립트 사용

---

## 추가로 필요한 정보 (분석 전 결정)

### 1. Target (Vestibule) Residue 범위

**필요한 정보:**
- GLUT1 vestibule 입구를 둘러싼 잔기 번호
- 예: `resid 50-80 150-190 280-320`

**확인 방법:**
```bash
# 00_inspect_selections.py 실행 후 protein residue 범위 확인
python 00_inspect_selections.py step3_input.psf
```

**참고 자료:**
- GLUT1 구조 논문
- 기존 docking/MD 연구 결과
- Vestibule 정의 (보통 extracellular gate 주변)

---

### 2. Tripod의 3개 Glucose Resid

**필요한 정보:**
- 말단 Glucose 3개의 residue ID
- 예: `resid 900, 901, 902`

**확인 방법:**
```bash
# 00_inspect_selections.py 실행 후 sugar residue 찾기
python 00_inspect_selections.py step3_input.psf
```

**출력 예시:**
```
[Residues that look like sugars]
  resid 900    resname GLC      natoms 24
  resid 901    resname GLC      natoms 24
  resid 902    resname GLC      natoms 24
```

---

## 파일 크기 예상

### Topology (.psf)
- **크기**: 수 MB ~ 수십 MB
- **Tripod 시스템**: ~10-50 MB (수용액 포함)

### Trajectory (.dcd)
- **크기**: 프레임당 ~1-5 MB
- **20ns (50,000 프레임)**: ~50-250 GB
- **압축 가능**: `.xtc` 형식으로 변환 시 ~10-50 GB

---

## 체크리스트

분석 시작 전 확인:

- [ ] Topology 파일 존재 확인
- [ ] Trajectory 파일 존재 확인
- [ ] 파일 크기 정상 (DCD > 1 GB)
- [ ] MDAnalysis 설치 확인 (`conda activate Drug-MD`)
- [ ] Target residue 범위 결정
- [ ] Glucose resid 3개 확인

---

## 다음 단계

### 1. 파일 확인
```bash
cd ~/projects/Drug/2026-01-22_Glycogate_Tripod
find . -name "*.psf" -o -name "*.dcd" | head -20
```

### 2. 선택자 확인
```bash
cd ~/projects/Drug/2026-01-23_AnchorR/scripts
python 00_inspect_selections.py \
    ~/projects/Drug/2026-01-22_Glycogate_Tripod/data/solution_builder/openmm/step3_input.psf
```

### 3. 분석 실행
```bash
python analyze_tripod_anchor_reach.py \
    ~/projects/Drug/2026-01-22_Glycogate_Tripod/data/solution_builder/openmm/step3_input.psf \
    ~/projects/Drug/2026-01-22_Glycogate_Tripod/results/md_3arm_tripod/step5_production.dcd
```

---

## 문제 해결

### "파일을 찾을 수 없습니다"
```bash
# 전체 프로젝트에서 PSF 파일 찾기
find ~/projects/Drug/2026-01-22_Glycogate_Tripod -name "*.psf"

# 전체 프로젝트에서 DCD 파일 찾기
find ~/projects/Drug/2026-01-22_Glycogate_Tripod -name "*.dcd"
```

### "DCD 파일이 너무 작음"
- MD 시뮬레이션이 완료되지 않았을 수 있음
- 진행 상황 확인 필요

### "PSF 파일이 없음"
- Solution Builder 단계가 완료되지 않았을 수 있음
- `step3_pbcsetup.oldpsf`를 `step3_input.psf`로 복사 가능
