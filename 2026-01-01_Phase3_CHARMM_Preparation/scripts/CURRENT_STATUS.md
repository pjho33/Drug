# Phase 3 시뮬레이션 현재 상태

## ❌ 발생한 문제들

### 1. PDB 포맷 제한
- OpenMM PDBFile은 원자 번호 99,999까지만 지원
- 현재 시스템: 113,807 atoms (실험군), 110,392 atoms (대조군)
- **해결 불가**: PDB 표준 포맷 제한

### 2. Parameter 로딩 문제
- HGTIP3 atom type 누락 (대조군)
- toppar.str include 경로 문제

### 3. PSF-PDB 불일치
- PSF에 Tripod topology 없음
- PDB에만 Tripod 좌표 있음

## 🎯 근본 원인

**기존 시스템에 Tripod를 "끼워넣기"하는 방식의 한계**

CHARMM-GUI는 처음부터 모든 구성요소를 포함하여 시스템을 생성하도록 설계됨

## ✅ 해결 방법

### 옵션 A: CHARMM-GUI Ligand Reader 사용 (권장)

**장점**:
- 완전히 일관된 시스템 생성
- PSF, PDB, toppar 모두 자동 생성
- 검증된 방식

**단점**:
- 2-3시간 대기 필요
- 재제출 필요

**진행 방법**:
1. CHARMM-GUI Ligand Reader 접속
2. Tripod SDF + GLUT1 PDB 업로드
3. Membrane Builder 설정
4. 다운로드 후 즉시 실행

### 옵션 B: 원본 CHARMM-GUI 스크립트 직접 사용

기존 CHARMM-GUI의 `openmm_run.py`를 그대로 사용 (Tripod 없이)

**장점**:
- 즉시 실행 가능
- 검증된 시스템

**단점**:
- Tripod 없음 (원래 목표 달성 불가)

## 💡 추천

**CHARMM-GUI Ligand Reader로 재생성**

이미 복합체 PDB가 있으므로:
- `glut1_tripod_complex.pdb` (실험군)
- `glut1_tripod_complex_control.pdb` (대조군)

이것들을 CHARMM-GUI Membrane Builder에 제출하면 됩니다.

## 📝 다음 단계

어떻게 진행할까요?

1. **CHARMM-GUI 재제출** (2-3시간 대기, 확실함)
2. **Tripod 없이 원본 시스템 실행** (즉시, 하지만 목표 미달성)
3. **다른 방법 시도** (예: GROMACS, NAMD 등)
