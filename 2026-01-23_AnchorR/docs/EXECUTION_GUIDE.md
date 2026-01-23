# AnchorR 분석 실행 가이드

## 준비 사항

### 필수 파일
1. **Topology 파일**: `.psf` (CHARMM-GUI/OpenMM)
2. **Trajectory 파일**: `.dcd`, `.xtc`, `.nc`

예시:
- `step5_production.psf`
- `step5_production.dcd` 또는 `production.dcd`

### 필수 소프트웨어
```bash
conda activate Drug-MD
# MDAnalysis, numpy 필요
```

## 실행 순서

### Step 1: 선택자 확인

```bash
cd ~/projects/Drug/2026-01-23_AnchorR/scripts
python 00_inspect_selections.py /path/to/step5_production.psf /path/to/production.dcd
```

**출력에서 확인할 것:**
1. **Glucose resid 3개** - Tripod의 3개 말단
   - 예: resid 900, 901, 902
   - resname이 GLC, BGLC, AGLC 등으로 나타날 수 있음

2. **Protein residue 범위** - Target (vestibule) 영역
   - 예: 50-80, 150-190, 280-320
   - GLUT1 vestibule 주변 잔기

### Step 2: 파라미터 설정

`analyze_tripod_anchor_reach.py` 파일을 열어서 수정:

```python
# 1. Target region (Step 1에서 확인한 vestibule residue 범위)
TARGET_SEL = "protein and resid 50-80 150-190 280-320"

# 2. Tripod 3 arms (Step 1에서 확인한 Glucose resid 3개)
ARM_SELS = [
    "resname GLC and resid 900",  # arm A
    "resname GLC and resid 901",  # arm B
    "resname GLC and resid 902",  # arm C
]

# 3. Cutoffs (필요시 조정)
ANCHOR_CUTOFF = 3.5  # Å
REACH_CUTOFFS = [3.5, 5.0, 8.0]  # Å
```

### Step 3: 본 분석 실행

```bash
python analyze_tripod_anchor_reach.py /path/to/step5_production.psf /path/to/production.dcd
```

**예상 실행 시간:**
- 10,000 프레임: ~5-10분
- 50,000 프레임: ~30-60분

## 결과 해석

### 출력 예시

```
=== Anchor Summary ===
Anchor frames (best arm within 3.5Å): 4523 / 10000 = 0.452

=== Conditional Reach Stats (given Anchor=True) ===
[rc= 3.5Å] P(>=2 arms | Anchor) = 0.123  |  P(3 arms | Anchor) = 0.008
[rc= 5.0Å] P(>=2 arms | Anchor) = 0.234  |  P(3 arms | Anchor) = 0.045
[rc= 8.0Å] P(>=2 arms | Anchor) = 0.456  |  P(3 arms | Anchor) = 0.123

=== Dwell Time ===
[rc= 5.0Å] events=234  mean=12.3 ps  p95=45.6 ps
```

### 해석 기준

| 지표 | 값 | 의미 |
|------|-----|------|
| Anchor 빈도 | > 30% | 기본 접근 가능 ✅ |
| P(>=2 arms \| Anchor) | > 10% | **Multivalency 가능** ⭐ |
| P(>=2 arms \| Anchor) | > 20% | **매우 우수** ⭐⭐ |
| P(3 arms \| Anchor) | 0-5% | 정상 (geometry 제약) |
| Dwell time | > 10 ps | 안정적 결합 가능 |

### 성공 기준

**PEG24 Tripod 설계 성공 판단:**
- ✅ P(>=2 arms | Anchor) > 10% at 5.0Å cutoff
- ✅ Dwell time > 10 ps
- ✅ Anchor 빈도 > 30%

## 출력 파일

### `tripod_anchor_reach.npz`

저장된 데이터:
- `arm_to_target`: (n_frames, 3) - 각 팔의 타겟까지 거리
- `anchor_arm`: (n_frames,) - 각 프레임의 anchor arm index
- `anchor_dist`: (n_frames,) - anchor arm의 거리
- `anchor_series`: (n_frames,) - anchor 여부 (True/False)
- `reach_cutoffs`: REACH_CUTOFFS 배열

### 추가 분석 (선택)

```python
import numpy as np
import matplotlib.pyplot as plt

data = np.load('tripod_anchor_reach.npz')
arm_to_target = data['arm_to_target']
anchor_series = data['anchor_series'].astype(bool)

# Anchor 프레임에서의 거리 분포
anchor_dists = arm_to_target[anchor_series]

plt.figure(figsize=(10, 6))
plt.hist(anchor_dists[:, 0], bins=50, alpha=0.5, label='Arm 0')
plt.hist(anchor_dists[:, 1], bins=50, alpha=0.5, label='Arm 1')
plt.hist(anchor_dists[:, 2], bins=50, alpha=0.5, label='Arm 2')
plt.xlabel('Distance to Target (Å)')
plt.ylabel('Count')
plt.legend()
plt.savefig('anchor_distance_distribution.png')
```

## 문제 해결

### 오류: "TARGET_SEL matched 0 atoms"
- Step 1의 `00_inspect_selections.py` 출력을 다시 확인
- Protein residue 범위가 올바른지 확인

### 오류: "ARM_SELS[i] matched 0 atoms"
- Glucose resid가 올바른지 확인
- resname이 GLC가 아닐 수 있음 (BGLC, AGLC 등)

### 실행 시간이 너무 길 경우
- Trajectory를 일부만 사용 (예: 처음 10,000 프레임)
- `u.trajectory[::10]`로 10 프레임마다 샘플링

## 다음 단계

1. ✅ 분석 완료 후 결과 해석
2. PEG16 vs PEG24 비교 분석
3. 다양한 REACH_CUTOFF 값으로 sensitivity 분석
4. Dwell time 분포 상세 분석
5. 시각화 및 논문 figure 작성
