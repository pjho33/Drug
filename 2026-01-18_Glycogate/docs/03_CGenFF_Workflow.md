# CGenFF 파라미터 생성 워크플로우

**작성일:** 2026-01-18  
**목적:** 1-Arm TRIS-PEG24-β‑L‑Glc CGenFF 파라미터 생성

---

## 🎯 개요

CGenFF (CHARMM General Force Field)를 사용하여 1-arm 리간드의 force field 파라미터를 생성합니다.

**입력:** MOL2 또는 PDB 파일  
**출력:** .str 파일 (CHARMM stream file)

---

## 1. 구조 파일 준비

### 1.1 SMILES → 3D 구조

```bash
cd /home/pjho3/projects/Drug/2026-01-18_Glycogate
python scripts/01_validate_and_generate_3d.py
```

**생성 파일:**
- `data/1arm_peg24_glc.mol2` - CGenFF 입력용
- `data/1arm_peg24_glc.pdb` - 시각화용
- `data/1arm_peg24_glc.sdf` - 백업

### 1.2 구조 확인

```bash
# PyMOL로 시각화
pymol data/1arm_peg24_glc.pdb

# 또는 VMD
vmd data/1arm_peg24_glc.pdb
```

**확인 사항:**
- [ ] PEG24 체인이 올바르게 연결되었는지
- [ ] Glucose 입체화학 (β-L)
- [ ] TRIS 중심 구조
- [ ] 원자 간 거리 정상 (결합 길이)

---

## 2. CGenFF 파라미터 생성

### 옵션 1: CGenFF Web Server (권장)

**URL:** https://cgenff.umaryland.edu/

**절차:**
1. 계정 생성/로그인
2. MOL2 파일 업로드 (`1arm_peg24_glc.mol2`)
3. Submit
4. 결과 다운로드 (.str 파일)

**장점:**
- 간편함
- 최신 버전 사용
- 자동 penalty score 계산

**단점:**
- 파일 크기 제한 (보통 100 원자 이하)
- 1-arm (~60 원자)는 가능, 3-arm은 불가능

---

### 옵션 2: 로컬 CGenFF 프로그램

**필요:**
- CGenFF 프로그램 (라이센스 필요)
- CHARMM36 파라미터 파일

**명령:**
```bash
# CGenFF 실행 (예시)
cgenff data/1arm_peg24_glc.mol2 -o data/1arm_peg24_glc.str
```

---

### 옵션 3: ParamChem (대안)

**URL:** https://www.paramchem.org/

CGenFF와 유사한 웹 서비스

---

## 3. Penalty Score 분석

### 3.1 .str 파일 확인

```bash
# Penalty score 확인
grep "PENALTY" data/1arm_peg24_glc.str | sort -k2 -n -r | head -20
```

**출력 예시:**
```
! PENALTY  45.00  dihedral  C1-C2-O3-C4
! PENALTY  12.50  angle     C1-O2-C3
! PENALTY   8.00  bond      C1-C2
```

### 3.2 Penalty Score 기준

| Score | 평가 | 조치 |
|-------|------|------|
| < 10 | 우수 | 사용 가능 |
| 10-50 | 양호 | 확인 필요 |
| 50-100 | 주의 | 검증 필수 |
| > 100 | 위험 | 재파라미터화 권장 |

### 3.3 주요 확인 항목

**High Penalty Dihedral:**
- PEG backbone (C-C-O-C)
- Triazole 연결부
- Amide 연결부
- O-glycosidic bond

**대응:**
1. Score < 50: 그대로 사용
2. Score 50-100: 1-arm MD로 검증
3. Score > 100: FFParam으로 재최적화

---

## 4. 파라미터 파일 구조

### .str 파일 구성

```
* CGenFF stream file for 1arm_peg24_glc
*

! Atom types
MASS   -1  CG331   12.01100  ! aliphatic C
MASS   -1  OG301   15.99940  ! ether oxygen
...

! Residue topology
RESI 1ARM        0.000  ! 1-arm ligand
GROUP
ATOM C1   CG331  -0.270
ATOM H1   HGA3    0.090
...

! Bonds
BOND C1  C2
BOND C2  O1
...

! Angles
ANGLE C1  C2  O1
...

! Dihedrals
DIHE C1  C2  O1  C3
! PENALTY  45.00
...

! Impropers
IMPR ...

! Non-bonded parameters
NONBONDED
CG331  0.0  -0.0780  2.050
...

END
```

---

## 5. CHARMM-GUI 통합

### 5.1 리간드 업로드

**CHARMM-GUI Solution Builder:**
1. Upload PDB/MOL2
2. Upload .str 파일
3. Topology/Parameter 자동 적용

### 5.2 시스템 구축

```
1-arm 리간드
    +
TIP3P 물
    +
150 mM NaCl
```

**Box size:** 12 × 12 × 12 nm³

---

## 6. 검증 체크리스트

### 파라미터 생성 후

- [ ] .str 파일 생성 확인
- [ ] Penalty score < 100 (전체)
- [ ] High penalty dihedral 확인 (> 50)
- [ ] Atom type 할당 확인
- [ ] Charge 중성 확인 (NH2 있으므로 +1 또는 중성)

### CHARMM-GUI 후

- [ ] Topology 오류 없음
- [ ] 에너지 최소화 성공
- [ ] 초기 구조 정상

---

## 7. 문제 해결

### 문제 1: 파일 크기 초과 (Web Server)

**증상:** "File too large" 오류

**해결:**
1. 분할 접근
   - PEG24만 따로
   - Glucose만 따로
   - 수동 병합
2. 로컬 CGenFF 사용
3. ParamChem 시도

---

### 문제 2: High Penalty Score (> 100)

**증상:** 특정 dihedral penalty > 100

**해결:**
1. **FFParam 사용**
   ```bash
   # QM 계산으로 dihedral 재파라미터화
   # Gaussian/ORCA 필요
   ```

2. **Force Field Toolkit (VMD)**
   - GUI 기반 파라미터 최적화
   - QM 데이터 fitting

3. **수동 조정**
   - .str 파일 직접 수정 (비권장)

---

### 문제 3: Glucose 입체화학 오류

**증상:** β-L 대신 β-D로 인식

**해결:**
1. SMILES 재확인
2. 3D 구조 수동 수정 (Avogadro/ChemDraw)
3. CHARMM36 당 파라미터 직접 사용

---

## 8. 다음 단계

### CGenFF 성공 시

1. **CHARMM-GUI로 시스템 구축**
   - Solution Builder
   - 1-arm + 물 + 이온

2. **OpenMM으로 변환**
   - ParmEd 사용
   - 또는 CHARMM-GUI OpenMM 출력

3. **1-arm MD 시뮬레이션**
   - 200-500 ns
   - 3 replica

### 재파라미터화 필요 시

1. **High penalty dihedral 목록 작성**
2. **QM 계산 준비**
   - Gaussian input 생성
   - Dihedral scan
3. **FFParam으로 fitting**
4. **재검증**

---

## 9. 참고 자료

### CGenFF 문서
- CGenFF paper: Vanommeslaeghe et al., J. Comput. Chem. 2010
- CHARMM36 documentation

### 도구
- **CGenFF Server:** https://cgenff.umaryland.edu/
- **ParamChem:** https://www.paramchem.org/
- **FFParam:** https://github.com/Acellera/ffparam
- **VMD Force Field Toolkit:** https://www.ks.uiuc.edu/Research/vmd/plugins/fftk/

### 예제
- PEG force field benchmarks
- Carbohydrate CGenFF examples

---

**작성자:** Cascade AI  
**최종 수정:** 2026-01-18
