# Phase 3 시뮬레이션 준비 상태

## ✅ 완료된 작업

1. **PDB 병합 완료**
   - 실험군: `experimental/step5_input_with_tripod.pdb` (113,807 atoms)
   - 대조군: `control/step5_input_with_tripod.pdb` (110,392 atoms)
   - Tripod 48 atoms 추가됨

2. **파일 준비**
   - Tripod topology: `trp.rtf`
   - Tripod parameters: `trp.prm`

---

## 🔧 현재 문제: PSF 업데이트

**문제**: 기존 PSF 파일에 Tripod topology를 추가해야 함

**해결 방법 3가지**:

### 옵션 1: CHARMM-GUI 원본 스크립트 활용 (권장)
기존 CHARMM-GUI의 `openmm_run.py`를 수정하여 사용
- Tripod 포함 PDB 사용
- `toppar.str`에 `trp.rtf`, `trp.prm` 추가
- PSF 없이 PDB만으로 시뮬레이션 가능

**장점**: 
- CHARMM-GUI가 검증한 방식
- PSF 수동 편집 불필요
- 즉시 실행 가능

**단점**:
- CHARMM-GUI 스크립트 구조 이해 필요

---

### 옵션 2: OpenMM Modeller 사용
OpenMM의 Modeller로 topology 자동 생성
- PSF 없이 PDB + ForceField로 시스템 생성
- Tripod residue template 등록 필요

**장점**:
- 유연함
- Python 코드로 완전 제어

**단점**:
- Tripod template 수동 등록 복잡
- CHARMM force field 호환성 확인 필요

---

### 옵션 3: PSF 수동 편집
기존 PSF에 Tripod 원자/결합 정보 직접 추가

**장점**:
- 기존 시스템 완전 유지

**단점**:
- 매우 복잡하고 오류 가능성 높음
- 권장하지 않음

---

## 💡 추천: 옵션 1 (CHARMM-GUI 스크립트 활용)

### 작업 순서

1. **toppar.str 수정**
```bash
# 기존 toppar.str 끝에 추가
echo "../trp.rtf" >> toppar.str
echo "../trp.prm" >> toppar.str
```

2. **openmm_run.py 수정**
```python
# PDB 경로를 Tripod 포함 버전으로 변경
-c step5_input.pdb
→ -c step5_input_with_tripod.pdb
```

3. **실행**
```bash
cd experimental/
python ../openmm_run.py -i step6.1_equilibration.inp \
    -p step5_input.psf \
    -c step5_input_with_tripod.pdb \
    -t toppar.str
```

---

## 🚀 다음 단계

어떤 방법으로 진행할까요?

1. **옵션 1** - CHARMM-GUI 스크립트 수정 (권장, 10분)
2. **옵션 2** - OpenMM Modeller (유연, 30분)
3. **옵션 3** - PSF 수동 편집 (복잡, 비권장)

선택하시면 해당 방법으로 즉시 진행하겠습니다.
