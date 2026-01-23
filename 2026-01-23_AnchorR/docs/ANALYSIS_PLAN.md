# AnchorR 분석 계획: Tripod Multivalency with Anchor Constraint

## 연구 목적

**TRIS-(PEG24-L-Glucose)₃ 구조에서 한 다리가 GLUT1에 고정(anchor)되었을 때, 나머지 다리들의 도달 거리 및 동시 결합(multivalency) 가능성 분석**

## 핵심 질문

> **조건부 분포 분석**: 다리 1개가 고정되어 있을 때, 나머지 2개 다리가 얼마나/얼마나 자주 타겟에 도달하는가?

### 주요 지표

1. **CRP1 (Conditional Reach Probability)**: P(≥1개 다리 추가 도달 | Anchor)
2. **CMP2 (Conditional Multivalency Probability)**: P(≥2개 다리 동시 도달 | Anchor) ← **핵심 지표**
3. **CMP3**: P(3개 모두 도달 | Anchor) - 참고용

## 분석 시나리오 3단계

### (A) Anchor 정의

**Anchor 조건**: 다리 A의 Glucose 말단이 타겟 residue set과 **최소 거리 ≤ 3.5 Å**

```
Anchor(frame) = True / False
```

### (B) Anchor 상태에서 나머지 다리의 거리 측정

Anchor가 True인 프레임만 필터링:
- 다리 B, C 각각의 end-to-end distance
- 또는 target residue까지의 min distance

### (C) 동시성(Multivalency) 판단

| 지표 | 의미 |
|------|------|
| ≥1 non-anchor arm reach | "두 번째 다리라도 뻗나?" |
| ≥2 non-anchor arms reach | "진짜 multivalency 되나?" ← **핵심** |
| ≥3 arms reach | "거의 이벤트 수준" |

## 핵심 지표 정의

### 🎯 지표 1: Conditional Reach Probability (CRP1)

```
CRP1 = (# frames where Anchor=True AND ≥1 other arm reaches)
       / (# frames where Anchor=True)
```

**의미**: 한 개 붙으면, 두 번째도 따라올 가능성

### 🎯 지표 2: Conditional Multivalency Probability (CMP2) ⭐

```
CMP2 = (# frames where Anchor=True AND ≥2 other arms reach)
       / (# frames where Anchor=True)
```

**의미**: **설계의 성패를 판단하는 핵심 지표**

### 🎯 보조 지표: Dwell Time

Anchor 상태에서 CMP2 상태가 **연속으로 유지되는 평균 시간**
- 짧은 깜빡임 vs 안정 결합 구분

## 파라미터 설정

### 필수 파라미터

```python
# 1. Target region (vestibule 주변 residue set)
TARGET_SEL = "protein and resid 50-80 150-190 280-320"

# 2. Tripod 3 arms (말단 glucose selection)
ARM_SELS = [
    "resname GLC and resid 900",  # arm A
    "resname GLC and resid 901",  # arm B
    "resname GLC and resid 902",  # arm C
]

# 3. Cutoffs (Å)
ANCHOR_CUTOFF = 3.5
REACH_CUTOFFS = [3.5, 5.0, 8.0]  # 여러 임계값으로 스윕
```

## 해석 가이드

| 결과 | 해석 |
|------|------|
| Anchor 빈도 높음 | 기본 접근 가능 |
| CRP1 ↑ | 두 번째 다리 가능성 있음 |
| **CMP2 > 5–10%** | **설계 성공 가능성 매우 높음** ⭐ |
| CMP3 ≈ 0 | 정상 (3개 동시는 geometry가 거의 완벽해야 함) |
| dwell time > 수 ns | 실제 결합 가능 |

## 중요 노트

### ✅ 올바른 접근
- **주 지표**: CMP2 (2개 이상 동시 도달)
- **보조 지표**: CMP3 (3개 동시)
- CMP3가 0이어도 실패가 아님

### ❌ 잘못된 접근
- 3개 동시를 목표로 판단 → 과도한 실패 판정
- Anchor 조건 없이 자유 상태 평균 계산 → 의미 없음

## 워크플로우

1. **선택자 확인**: `00_inspect_selections.py` 실행
   - Glucose resid 3개 찾기
   - Target residue 범위 확인

2. **본 분석**: `analyze_tripod_anchor_reach.py` 실행
   - Anchor 조건부 도달 거리 계산
   - CRP1, CMP2, CMP3 계산
   - Dwell time 분석

3. **결과 해석**
   - CMP2 > 5-10%: 성공 가능성 높음
   - Dwell time 분포 확인
   - PEG 길이별 비교 (PEG16 vs PEG24)

## 예상 결과

### PEG24 (현재)
- CMP2: 10-20% 예상
- CMP3: 0-2% (정상)
- Dwell time: 수 ns

### PEG16 (비교)
- CMP2: 5-10% (더 낮음)
- CMP3: ~0%

## 다음 단계

1. Tripod MD 시뮬레이션 완료 확인
2. 선택자 확인 스크립트 실행
3. 본 분석 스크립트 실행
4. 결과 시각화 및 해석
5. PEG 길이별 비교 분석
