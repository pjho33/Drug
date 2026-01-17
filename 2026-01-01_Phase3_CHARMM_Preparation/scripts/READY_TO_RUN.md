# Phase 3 시뮬레이션 준비 완료

## ✅ 확인 완료

1. **PDB 병합 완료**
   - 실험군: 113,807 atoms (GLUT1 + Glycans + Tripod + membrane)
   - 대조군: 110,392 atoms (GLUT1 + Tripod + membrane)

2. **Tripod 파일 확인**
   - Tripod.sdf: 복합체와 동일 (MD5 일치)
   - 48 heavy atoms (수소는 CHARMM이 자동 추가)
   - PEG6 버전 사용

3. **Topology/Parameters**
   - trp.rtf: Tripod topology
   - trp.prm: Tripod parameters

---

## 🎯 시뮬레이션 실행 방법

### 옵션 1: CHARMM-GUI 원본 스크립트 활용 (권장)

**장점**: 검증된 방식, PSF 수동 편집 불필요

**단계**:

1. **toppar.str 수정**
```bash
# 실험군
cd /home/pjho3tr/Downloads/charmm-gui-6750265216membranebuilder/openmm
cp toppar.str toppar.str.backup
echo "../../projects/Drug/phase3_with_tripod/trp.rtf" >> toppar.str
echo "../../projects/Drug/phase3_with_tripod/trp.prm" >> toppar.str

# 대조군
cd /home/pjho3tr/Downloads/charmm-gui-6704990786대조군/openmm
cp toppar.str toppar.str.backup
echo "../../projects/Drug/phase3_with_tripod/trp.rtf" >> toppar.str
echo "../../projects/Drug/phase3_with_tripod/trp.prm" >> toppar.str
```

2. **Tripod 포함 PDB 복사**
```bash
# 실험군
cp /home/pjho3tr/projects/Drug/phase3_with_tripod/experimental/step5_input_with_tripod.pdb \
   /home/pjho3tr/Downloads/charmm-gui-6750265216membranebuilder/openmm/

# 대조군
cp /home/pjho3tr/projects/Drug/phase3_with_tripod/control/step5_input_with_tripod.pdb \
   /home/pjho3tr/Downloads/charmm-gui-6704990786대조군/openmm/
```

3. **실행 스크립트 작성**
```bash
# phase3_with_tripod/run_experimental_gpu0.py
# phase3_with_tripod/run_control_gpu1.py
```

---

### 옵션 2: OpenMM Modeller 사용

PSF 없이 PDB + ForceField로 시스템 생성

**장점**: 유연함
**단점**: Tripod template 등록 필요

---

## 📝 다음 작업

어떤 방법으로 진행할까요?

1. **옵션 1 실행** - CHARMM-GUI 스크립트 수정 및 실행
2. **옵션 2 실행** - OpenMM Modeller로 새로 생성

선택하시면 즉시 준비하겠습니다.
