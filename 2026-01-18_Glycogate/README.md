# Glycogate: Tripod-PEG24-Glucose 시스템 검증

**시작일:** 2026-01-18

## 📋 프로젝트 개요

**목적:** Tripod-PEG24-Triazole-Amide-Glucose 시스템의 Force Field 검증 및 물성 확인

**시스템 구성:**
- TRIS 분기점 (Tripod core)
- PEG24 링커 (-(CH2-CH2-O)24-)
- Triazole 연결 (1,2,3-triazole)
- Amide 연결 (-CO-NH-)
- Glucose (β-D-Glucose)

**핵심 질문:**
1. PEG24가 물에서 얼마나 펼쳐지는가?
2. CHARMM36 + CGenFF 조합이 신뢰할 만한가?
3. RBC 회피 vs 종양내피 접근에 필요한 r_eff, r_tail은?

---

## 📁 폴더 구조

```
2026-01-18_Glycogate/
├── scripts/      # MD 시뮬레이션 및 분석 스크립트
├── data/         # 구조 파일, 파라미터 파일
├── results/      # MD trajectory, 분석 결과
├── docs/         # 전략 문서, 프로토콜
└── README.md     # 이 파일
```

---

## 📚 문서

1. **Force Field 전략** (`docs/01_Force_Field_Strategy.md`)
   - CHARMM36 + CGenFF 선택 근거
   - 검증 체크리스트
   - 문제 해결 방법

2. **1-Arm 검증 프로토콜** (`docs/02_1-Arm_Validation_Protocol.md`)
   - 시스템 구축 방법
   - MD 시뮬레이션 설정
   - 분석 프로토콜
   - 판정 기준

---

## 🔧 스크립트 목록

| 파일명 | 설명 |
|--------|------|
| `run_example.py` | 예시 스크립트 (수정 필요) |

*(작업 진행 중)*

---

## 🚀 워크플로우

### Phase 1: 1-Arm 검증 (현재)

1. **구조 준비**
   - 1-arm 모델 구축 (TRIS-stub + PEG24 + Triazole + Amide + Glc)
   - CGenFF 파라미터 생성
   - Penalty score 확인

2. **MD 시뮬레이션**
   - 수용액 (TIP3P + 150 mM NaCl)
   - NPT, 300 K, 200-500 ns
   - 3 replica (다른 seed)

3. **분석**
   - End-to-end 거리 분포
   - Radius of Gyration
   - Tail 확률 (P(R > 6 nm), P(R > 8 nm))
   - Autocorrelation time
   - Dihedral 분포

4. **판정**
   - CGenFF 사용 가능 여부 결정
   - 3-arm 진행 또는 재최적화

### Phase 2: 3-Arm 시스템 (예정)

1. Tripod 전체 시스템 구축
2. 팔-팔 상호작용 분석
3. 최종 물성 확정

---

## 📊 예상 결과

### 정상 범위 (문헌 기준)

- **End-to-end 거리:** 2-4 nm (주), 6-9 nm (꼬리)
- **Rg:** 1.5-2.5 nm
- **P(R > 6 nm):** 0.05-0.20
- **P(R > 8 nm):** 0.001-0.05

### 활용

- RBC 회피 확률 계산
- 종양내피 접근 모델링
- 약동학 파라미터 도출

## 📊 데이터

- **입력:** `data/` 폴더
- **출력:** `results/` 폴더
- **대용량 파일:** Git에서 제외 (`.gitignore` 참조)

## 🚀 실행 방법

```bash
cd scripts/
python run_*.py
```

## 📝 노트

- 추가 정보 및 메모

---

**최종 수정:** 2026-01-18
