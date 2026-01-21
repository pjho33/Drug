# Glycogate Bipod: TRIS-(PEG24-L-Glucose)₂ 시스템 검증

**시작일:** 2026-01-20

## 📋 프로젝트 개요

**목적:** TRIS-(PEG24-L-Glucose)₂ Bipod 시스템의 Force Field 검증 및 물성 확인

**시스템 구성:**
- TRIS 분기점 (Bipod core - 2개 팔)
- PEG24 링커 × 2 (-(CH2-CH2-O)24-)
- Triazole 연결 (1,2,3-triazole)
- Amide 연결 (-CO-NH-)
- L-Glucose × 2 (β-L-Glucose)

**1-Arm vs 2-Arm 비교:**
- **1-Arm (2026-01-18):** TRIS-(OH)₂-PEG24-L-Glc (단일 팔)
- **2-Arm (현재):** TRIS-OH-(PEG24-L-Glc)₂ (두 개 팔)

**핵심 질문:**
1. 2개 팔 사이의 상호작용이 있는가?
2. 팔-팔 간섭이 PEG 확장에 영향을 주는가?
3. 1-arm 대비 effective radius가 어떻게 변하는가?
4. CHARMM36 + CGenFF가 2-arm 시스템에서도 신뢰할 만한가?

---

## 📁 폴더 구조

```
2026-01-20_Glycogate_Bipod/
├── scripts/      # MD 시뮬레이션 및 분석 스크립트
├── data/         # 구조 파일, 파라미터 파일
├── results/      # MD trajectory, 분석 결과
├── docs/         # 전략 문서, 프로토콜
└── README.md     # 이 파일
```

---

## 🚀 워크플로우

### Phase 1: 2-Arm 시스템 구축

1. **구조 준비**
   - 2-arm 모델 구축 (TRIS-OH + 2×(PEG24 + Triazole + Amide + L-Glc))
   - CGenFF 파라미터 생성
   - Penalty score 확인

2. **MD 시뮬레이션**
   - 수용액 (TIP3P + 150 mM NaCl)
   - NPT, 300 K, 200-500 ns
   - 3 replica (다른 seed)

3. **분석**
   - **각 팔의 End-to-end 거리 분포**
   - **팔-팔 간 거리 및 각도**
   - Radius of Gyration (전체 시스템)
   - Tail 확률 (P(R > 6 nm), P(R > 8 nm))
   - Autocorrelation time
   - Dihedral 분포

4. **1-Arm과 비교**
   - PEG 확장도 차이
   - 팔 간 상호작용 정량화
   - Effective size 변화

---

## 📊 예상 결과

### 2-Arm 특성 (예상)

- **각 팔의 End-to-end 거리:** 1-arm과 유사하거나 약간 감소 (간섭 효과)
- **팔-팔 간 각도:** 90-180° (자유 회전)
- **전체 Rg:** 1-arm의 1.2-1.5배
- **팔 간 상호작용:** 약한 배제 부피 효과

### 1-Arm 대비 차이점

- **Effective radius:** 증가 (2개 팔)
- **RBC 회피:** 더 유리 (크기 증가)
- **종양내피 접근:** 약간 불리할 수 있음 (크기 증가)

---

## 📊 데이터

- **입력:** `data/` 폴더
- **출력:** `results/` 폴더
- **대용량 파일:** Git에서 제외 (`.gitignore` 참조)

---

## 🚀 실행 방법

```bash
# 환경 활성화
conda activate Drug-MD

# 스크립트 실행
cd scripts/
python run_*.py
```

---

## 📝 비교 분석

### 1-Arm (2026-01-18) vs 2-Arm (현재)

| 특성 | 1-Arm | 2-Arm (예상) |
|------|-------|--------------|
| 팔 개수 | 1 | 2 |
| Rg | 1.5-2.5 nm | 2.0-3.5 nm |
| 팔 간 상호작용 | 없음 | 약함-중간 |
| 시뮬레이션 비용 | 낮음 | 중간 |

---

## 🎯 다음 단계

1. 2-Arm 시스템 검증 완료
2. 1-Arm vs 2-Arm 비교 분석
3. 3-Arm (Tripod) 시스템으로 확장 결정

---

**최종 수정:** 2026-01-20
