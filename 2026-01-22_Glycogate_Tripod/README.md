# 2026-01-22 Glycogate Tripod Project

## 프로젝트 개요

**3-Arm (Tripod) Glycogate 시스템의 MD 시뮬레이션**

TRIS 중심에서 3개의 PEG24-L-Glucose 팔이 연결된 tripod 구조의 분자 동역학 시뮬레이션 프로젝트입니다.

### 시스템 구성
- **중심 구조:** TRIS (Tris(hydroxymethyl)aminomethane)
- **팔 구조:** 3 × (PEG24-L-Glucose)
- **총 팔 개수:** 3개 (Tripod)
- **용매:** 명시적 수용액 (Explicit water)
- **Force Field:** CHARMM36

## 폴더 구조

```
2026-01-22_Glycogate_Tripod/
├── README.md                    # 이 파일
├── data/                        # 입력 데이터 및 구조 파일
│   ├── 3arm_peg24_glc.smi      # 3-arm SMILES 구조
│   ├── Ligand reader/          # CHARMM-GUI Ligand Reader 파일
│   └── solution builder/       # CHARMM-GUI Solution Builder 파일
├── results/                     # MD 시뮬레이션 결과
│   └── md_3arm_tripod/         # Tripod MD 결과
├── scripts/                     # 실행 스크립트
│   ├── 09_run_openmm_simple.sh # MD 실행 스크립트
│   ├── start_md_background.sh  # 백그라운드 실행
│   └── check_md_progress.sh    # 진행 상황 확인
└── docs/                        # 문서 및 분석
```

## 워크플로우

### 1. SMILES 생성
- 3-arm tripod 구조의 SMILES 문자열 생성
- RDKit으로 수소 추가

### 2. CGenFF 파라미터 생성
- https://cgenff.umaryland.edu/
- SMILES 업로드 → `lig.rtf`, `lig.prm` 다운로드

### 3. Ligand Reader
- CHARMM-GUI Ligand Reader로 3D 구조 생성
- 에너지 최소화된 리간드 구조 생성

### 4. Solution Builder
- CHARMM-GUI Solution Builder로 수용액 시스템 구축
- 명시적 수용액 추가
- 이온 농도 설정

### 5. MD 시뮬레이션
- 에너지 최소화
- 평형화 (125 ps)
- Production MD (20 ns)
- CUDA 가속 사용

## 예상 결과

### 구조적 특성
- **Tripod 형태:** 3개 팔의 공간 배치
- **팔 간 상호작용:** 3개 팔 사이의 상호작용 패턴
- **회전 반경 (Rg):** Tripod의 전체 크기
- **SASA:** 용매 접근 가능 표면적

### 동역학적 특성
- **팔의 유연성:** 각 팔의 움직임
- **중심 안정성:** TRIS 중심의 안정성
- **글루코스 배향:** 3개 글루코스의 공간적 배치

### Bipod와 비교
- **크기:** Tripod > Bipod
- **안정성:** 3개 팔의 균형
- **용매 상호작용:** 더 큰 표면적

## 실행 방법

### 백그라운드 실행
```bash
cd ~/projects/Drug/2026-01-22_Glycogate_Tripod/scripts
bash start_md_background.sh
```

### 진행 상황 확인
```bash
bash check_md_progress.sh
```

### 로그 실시간 보기
```bash
tail -f ~/projects/Drug/2026-01-22_Glycogate_Tripod/scripts/md_run.log
```

## 시스템 요구사항

- **GPU:** NVIDIA CUDA 지원 (RTX 3090 권장)
- **메모리:** 16GB 이상
- **디스크:** 50GB 이상 여유 공간
- **소프트웨어:**
  - OpenMM 8.1.2
  - CUDA 13.1
  - Python 3.10
  - Conda 환경: Drug-MD

## 관련 프로젝트

- **1-Arm (Monopod):** `2026-01-18_Glycogate/`
- **2-Arm (Bipod):** `2026-01-20_Glycogate_Bipod/`
- **3-Arm (Tripod):** 이 프로젝트

## 다음 단계

1. 3-arm SMILES 생성
2. CGenFF 파라미터 생성
3. Ligand Reader 실행
4. Solution Builder 실행
5. 20ns MD 시뮬레이션
6. 결과 분석 및 1-arm, 2-arm과 비교

## 참고 사항

- Tripod 구조는 Bipod보다 복잡하여 시뮬레이션 시간이 더 길 수 있습니다
- CUDA 가속 필수 (CPU로는 매우 느림)
- 3개 팔의 대칭성을 유지하면서 구조 생성 필요
