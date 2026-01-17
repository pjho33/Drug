# GLUT1-Tripod 시뮬레이션 준비 완료

**작성일**: 2026-01-04  
**상태**: ✅ 시뮬레이션 준비 완료

---

## 🎯 목표

Glycocalyx 유무에 따른 Tripod 약물의 GLUT1 접근성 비교

---

## 📁 생성된 파일 (시뮬레이션용)

### 위치: `/home/pjho3/projects/Drug/scripts/`

### 1️⃣ **실험군 (Glycosylated GLUT1)**

| 파일 | 크기 | 원자 수 | 설명 |
|------|------|---------|------|
| **`glut1_tripod_complex.pdb`** | 587 KB | 7,702 | **실험군 복합체** (GLUT1 + Glycans + Tripod) |
| `glut1_glycosylated.pdb` | 584 KB | 7,654 | Glycosylated GLUT1만 |
| `tripod_aligned.pdb` | 3.9 KB | 48 | Tripod만 (정렬된 좌표) |

**특징**:
- Glycan 부착: Asn288, Asn88 (bi-antennary N-glycans with sialic acid)
- Tripod 위치: (584.18, -24.99, 202.84) Å
- 예상: Glycan이 Tripod 접근 방해

---

### 2️⃣ **대조군 (Non-Glycosylated GLUT1)**

| 파일 | 크기 | 원자 수 | 설명 |
|------|------|---------|------|
| **`glut1_tripod_complex_control.pdb`** | 543 KB | 7,114 | **대조군 복합체** (GLUT1 + Tripod, NO glycans) |
| `glut1_control.pdb` | 539 KB | 7,066 | Non-glycosylated GLUT1만 |

**특징**:
- Glycan 없음 (naked GLUT1)
- Tripod 위치: 실험군과 동일 (584.18, -24.99, 202.84) Å
- 예상: Tripod가 자유롭게 결합

---

### 3️⃣ **Topology & Parameters**

| 파일 | 크기 | 설명 |
|------|------|------|
| `trp.rtf` | 8.5 KB | Tripod topology (CHARMM format) |
| `trp.prm` | 8.3 KB | Tripod parameters (CGenFF) |

**Residue name**: `TRP`

---

## 🔬 시뮬레이션 설정 권장사항

### **공통 설정**

**1. Membrane Builder (CHARMM-GUI)**
- Input: 복합체 PDB 파일
- Membrane: POPC 또는 tumor cell composition
- Box size: 충분한 크기 (최소 120 Å × 120 Å)
- Water: TIP3P
- Ions: 0.15 M KCl

**2. MD Protocol**
```
Minimization: 5,000 steps
Equilibration: 
  - NVT: 100 ps (restraints on protein)
  - NPT: 500 ps (gradually release restraints)
Production: 100-200 ns
```

**3. Analysis**
- RMSD: Tripod position stability
- Distance: Tripod - glucose pocket
- Distance: Tripod - glycan (실험군만)
- Binding free energy: MM/PBSA or MM/GBSA
- Contact analysis: Tripod - GLUT1 residues

---

### **실험군 특이 설정**

**Glycan 고려사항**:
- Glycan flexibility: 높은 자유도
- Sialic acid: 음전하 (-1) → Tripod와 정전기적 반발
- 분석 추가: Glycan-Tripod 최소 거리

**예상 결과**:
- Tripod가 glycan에 의해 차단
- Binding affinity 감소
- Residence time 짧음

---

### **대조군 특이 설정**

**Naked GLUT1**:
- Glycan 없음 → 직접 접근 가능
- 표면 전하 다름 (sialic acid 없음)

**예상 결과**:
- Tripod가 glucose pocket에 안정적으로 결합
- Binding affinity 높음
- Residence time 길음

---

## 📊 비교 분석 계획

### **정량적 지표**

| 지표 | 실험군 (glycan) | 대조군 (no glycan) | 예상 차이 |
|------|----------------|-------------------|----------|
| **Binding ΔG** | 높음 (약한 결합) | 낮음 (강한 결합) | > 5 kcal/mol |
| **RMSD (Tripod)** | 높음 (불안정) | 낮음 (안정) | 2-3 Å |
| **Min distance (Tripod-Glycan)** | 5-10 Å | N/A | - |
| **Residence time** | 짧음 | 길음 | 10배 이상 |
| **Contact number** | 적음 | 많음 | 2배 이상 |

---

## 🎯 가설 검증

**가설**: Glycocalyx가 Tripod의 GLUT1 접근을 차단한다

**검증 방법**:
1. ✅ 실험군 vs 대조군 비교
2. ✅ 동일한 Tripod 초기 위치 사용
3. ✅ 동일한 시뮬레이션 조건

**성공 기준**:
- 실험군에서 Tripod-GLUT1 결합력 유의미하게 감소
- 대조군에서 Phase 2 수준의 결합력 유지
- Glycan-Tripod 거리가 접근 차단 증명

---

## 📝 다음 단계 (다른 컴퓨터에서)

### 1. **파일 전송**
```bash
# 필수 파일
glut1_tripod_complex.pdb          # 실험군
glut1_tripod_complex_control.pdb  # 대조군
trp.rtf                            # Topology
trp.prm                            # Parameters
```

### 2. **CHARMM-GUI Membrane Builder**
- 각 복합체를 별도로 업로드
- Membrane 추가
- Solvation & Ionization
- OpenMM 또는 GROMACS 입력 파일 생성

### 3. **시뮬레이션 실행**
```bash
# 실험군
cd experimental_group
python run_openmm.py

# 대조군
cd control_group
python run_openmm.py
```

### 4. **분석**
```bash
# RMSD, distance, binding energy 계산
python analyze_trajectory.py
```

---

## 🔍 검증 완료

### ✅ **구조 검증**

**실험군**:
- GLUT1: 7,654 atoms (glycans 포함)
- Tripod: 48 atoms
- Glycan sites: Asn288, Asn88
- 정렬 RMSD: < 0.01 Å

**대조군**:
- GLUT1: 7,066 atoms (glycans 없음)
- Tripod: 48 atoms (동일 위치)
- 정렬 RMSD: < 0.01 Å

### ✅ **좌표 검증**

**Tripod 중심 (두 군 동일)**:
- X: 584.18 Å
- Y: -24.99 Å
- Z: 202.84 Å

**Phase 2 MD 결과 기반**:
- 안정화된 결합 포즈
- Glucose pocket 내부 위치

---

## 📚 참고 자료

**Phase 2 결과**:
- `/media/pjho3/pjho backup/Drug 251228/Drug251231/projects251231Repliphase2/`

**CHARMM-GUI 원본**:
- 실험군: `/home/pjho3/다운로드/charmm-gui-6750265216membranebuilder/`
- 대조군: `/home/pjho3/다운로드/charmm-gui-6704990786대조군/`
- Ligand: `/home/pjho3/다운로드/charmm-gui-6753234611Ligand260104/`

---

## ✅ 체크리스트

- [x] Phase 2 Tripod 좌표 추출
- [x] CHARMM-GUI Ligand topology 확인
- [x] 실험군 GLUT1 구조 확인 (glycosylated)
- [x] 대조군 GLUT1 구조 확인 (non-glycosylated)
- [x] Kabsch alignment 수행
- [x] 실험군 복합체 생성
- [x] 대조군 복합체 생성
- [x] 파일 검증 완료
- [ ] CHARMM-GUI Membrane Builder 실행 (다른 컴퓨터)
- [ ] MD 시뮬레이션 실행
- [ ] 결과 분석

---

**준비 완료!** 🎉

모든 파일이 시뮬레이션 준비되었습니다.
다른 컴퓨터로 파일을 전송하여 시뮬레이션을 진행하세요.
