# CHARMM-GUI Ligand Reader 제출 가이드

## 목적
Tripod 리간드를 포함한 GLUT1 membrane 시스템을 생성하여 PSF/PDB/toppar 호환성 문제 해결

## 필요 파일

### 1. Protein Structure
- **Glycosylated GLUT1**: `/home/pjho3tr/Downloads/charmm-gui-6750265216membranebuilder/openmm/step5_input.pdb`
- **Control GLUT1**: `/home/pjho3tr/Downloads/charmm-gui-6704990786대조군/openmm/step5_input.pdb`

### 2. Ligand
- **Tripod SDF**: `tripod.sdf` (이 디렉토리에 복사됨)

### 3. Ligand 배치 좌표
- **Phase 2 결과**: `/home/pjho3tr/projects/Drug/results/phase2_rep1/prod_tripod_rep1_final.pdb`
- Tripod의 최종 결합 위치 참조용

## CHARMM-GUI 제출 절차

### Step 1: Ligand Reader 접속
https://www.charmm-gui.org/?doc=input/ligandrm

### Step 2: Ligand 업로드
1. "Upload Ligand" 섹션
2. `tripod.sdf` 파일 업로드
3. Ligand name: **TRIPOD**
4. Force field: **CGenFF** (CHARMM General Force Field)

### Step 3: Protein 업로드
1. "Upload PDB" 섹션
2. Glycosylated: `step5_input.pdb` (glycosylated system)
3. Control: `step5_input.pdb` (control system)
4. 각각 별도로 제출 (2개 시스템)

### Step 4: Ligand 위치 지정
**Option A: 자동 도킹**
- CHARMM-GUI의 AutoDock 사용
- Binding site: GLUT1 glucose binding pocket
- Grid center: 대략 (x, y, z) = (0, 0, 0) 근처 (PDB 좌표계)

**Option B: 수동 배치 (권장)**
1. Phase 2 결과에서 Tripod 좌표 추출
2. CHARMM-GUI에서 "Specify ligand position" 옵션
3. 좌표 직접 입력 또는 참조 PDB 업로드

### Step 5: Membrane Builder 설정
**이전 설정과 동일하게 유지**:
- Membrane type: **POPC** (또는 이전 사용한 지질)
- System size: **기존과 동일** (약 10 nm × 10 nm)
- Ion concentration: **0.15 M KCl**
- Water model: **TIP3P**

### Step 6: Glycosylation (Glycosylated system만)
- Asn45에 glycan 추가 (이전과 동일)
- Glycan type: 기존 설정 참조

### Step 7: OpenMM 옵션
- Output format: **OpenMM**
- Include: PSF, PDB, toppar files
- Simulation ready: Yes

## 예상 결과

각 시스템당 다음 파일들을 받게 됩니다:
```
charmm-gui-XXXXXXX/
├── openmm/
│   ├── step5_input.psf      # Tripod 포함 PSF
│   ├── step5_input.pdb      # Tripod 포함 PDB
│   ├── step5_input.crd      # 좌표 파일
│   ├── toppar.str           # Parameter 목록
│   ├── toppar/              # 모든 parameter 파일
│   ├── openmm_run.py        # 실행 스크립트
│   └── step6.*.inp          # 입력 파일들
```

## 다운로드 후 작업

1. **파일 압축 해제**
```bash
cd /home/pjho3tr/Downloads/
unzip charmm-gui-XXXXXXX-glyco-tripod.zip
unzip charmm-gui-YYYYYYY-control-tripod.zip
```

2. **MD 실행 스크립트 수정**
- GPU 설정
- 100 ns production run
- Checkpoint/resume 기능

3. **시뮬레이션 실행**
```bash
cd /home/pjho3tr/projects/Drug/phase3_glycosylation
bash run_both_gpus.sh
```

## 주의사항

1. **두 시스템을 별도로 제출**: Glycosylated와 Control을 각각 독립적으로
2. **Ligand 위치 일관성**: 두 시스템에서 Tripod 위치가 동일하도록
3. **Membrane 설정 일관성**: 지질 조성, 시스템 크기 동일하게
4. **Glycan 차이만 유지**: Glycosylated에만 Asn45 glycan

## 타임라인

- CHARMM-GUI 제출: **즉시**
- 처리 시간: **30분 ~ 2시간**
- 다운로드 후 MD 시작: **즉시**
- 100 ns MD 완료: **약 24-48시간** (GPU 성능에 따라)

## 문제 해결

만약 CHARMM-GUI에서 에러가 발생하면:
1. Ligand SDF 파일 형식 확인
2. Protein PDB에서 불필요한 HETATM 제거
3. Ligand 위치가 protein과 너무 가까운지 확인 (최소 2 Å 거리)
