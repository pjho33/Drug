# AnchorR Project (2026-01-23)

## 프로젝트 개요

**AnchorR 시스템 연구 프로젝트**

이 프로젝트는 AnchorR 시스템의 구조, 동역학 및 기능적 특성을 연구하기 위한 분자 동역학 시뮬레이션 및 분석 프로젝트입니다.

## 폴더 구조

```
AnchorR_2026-01-23/
├── README.md                    # 이 파일
├── data/                        # 입력 데이터 및 구조 파일
│   ├── structures/              # 분자 구조 파일
│   ├── parameters/              # Force field 파라미터
│   └── input_files/             # 시뮬레이션 입력 파일
├── results/                     # 시뮬레이션 결과
│   ├── md_trajectories/         # MD 궤적 파일
│   ├── analysis/                # 분석 결과
│   └── figures/                 # 그래프 및 시각화
├── scripts/                     # 실행 스크립트
│   ├── setup/                   # 시스템 준비 스크립트
│   ├── simulation/              # MD 실행 스크립트
│   └── analysis/                # 분석 스크립트
└── docs/                        # 문서 및 노트
    ├── notes/                   # 연구 노트
    ├── protocols/               # 실험 프로토콜
    └── references/              # 참고 문헌
```

## 시스템 요구사항

- **GPU:** NVIDIA CUDA 지원 (권장)
- **메모리:** 16GB 이상
- **디스크:** 충분한 여유 공간
- **소프트웨어:**
  - OpenMM 8.1.2
  - Python 3.10
  - Conda 환경: Drug-MD

## 워크플로우

### 1. 시스템 준비
- 구조 파일 준비
- Force field 파라미터 설정
- 시뮬레이션 입력 파일 생성

### 2. MD 시뮬레이션
- 에너지 최소화
- 평형화
- Production MD

### 3. 분석
- 구조적 분석
- 동역학적 분석
- 결과 시각화

## 실행 방법

### 시뮬레이션 실행
```bash
cd ~/projects/Drug/AnchorR_2026-01-23/scripts
# 실행 스크립트 추가 예정
```

### 분석 실행
```bash
cd ~/projects/Drug/AnchorR_2026-01-23/scripts/analysis
# 분석 스크립트 추가 예정
```

## 다음 단계

1. 시스템 구조 정의
2. Force field 파라미터 준비
3. MD 시뮬레이션 설정
4. 시뮬레이션 실행
5. 결과 분석

## 참고 사항

- 프로젝트 진행 상황은 `docs/notes/`에 기록
- 모든 스크립트는 재현 가능하도록 작성
- 결과 파일은 적절히 백업

## 관련 프로젝트

- Glycogate 시리즈 (1-arm, 2-arm, 3-arm)
- Phase2 Additional MD
- SDG Glycated GLUT1 Complex
