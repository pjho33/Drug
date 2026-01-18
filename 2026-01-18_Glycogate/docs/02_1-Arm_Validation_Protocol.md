# 1-Arm 검증 프로토콜

**작성일:** 2026-01-18  
**목적:** PEG24-Triazole-Amide-Glucose 1팔 시스템 검증

---

## 🎯 검증 목표

1. PEG24가 물에서 얼마나 펼쳐지는가?
2. CGenFF 파라미터가 신뢰할 만한가?
3. 3-arm 시스템 구축 전 물성 확정

---

## 1. 시스템 구축

### 1.1 분자 구조

```
[Cap/TRIS-stub]
    |
  PEG24 (-(CH2-CH2-O)24-)
    |
  Triazole (1,2,3-triazole)
    |
  Amide (-CO-NH-)
    |
  CH2
    |
  O-glycosidic bond
    |
  Glucose (β-D-Glucose)
```

### 1.2 구조 준비

**옵션 1: TRIS 분기점 포함**
- TRIS 중심 유지
- 나머지 2팔을 methyl로 캡핑

**옵션 2: 단순 캡 (권장)**
- TRIS 대신 작은 캡 (acetyl 등)
- 최소 시스템으로 검증

### 1.3 파라미터 생성

```bash
# CGenFF 파라미터 생성
# 1. SMILES/MOL2 파일 준비
# 2. CGenFF server 또는 로컬 도구 사용
# 3. .str 파일 생성
# 4. Penalty score 확인
```

**Penalty Score 기준:**
- < 10: 우수
- 10-50: 양호 (확인 필요)
- > 50: 재검토 필요 (특히 dihedral)

---

## 2. MD 시뮬레이션 설정

### 2.1 시스템 준비

| 항목 | 설정 |
|------|------|
| 용매 | TIP3P water |
| Box type | Cubic |
| Box size | 최소 1.5 nm 여유 |
| 이온 | 150 mM NaCl |
| 중화 | Yes |

**Box size 계산:**
```python
# 최대 end-to-end 거리 예상: ~9 nm
# Box size = 9 + 2*1.5 = 12 nm
# 권장: 12 × 12 × 12 nm³
```

### 2.2 시뮬레이션 파라미터

#### Energy Minimization
```
Algorithm: Steepest Descent
Max steps: 50,000
Force tolerance: 1000 kJ/mol/nm
```

#### NVT Equilibration
```
Duration: 100 ps
Temperature: 300 K
Thermostat: V-rescale (tau_t = 0.1 ps)
Constraints: h-bonds (LINCS)
Timestep: 2 fs
```

#### NPT Equilibration
```
Duration: 1 ns
Temperature: 300 K
Pressure: 1 bar
Barostat: Parrinello-Rahman (tau_p = 2.0 ps)
Thermostat: V-rescale (tau_t = 0.1 ps)
Constraints: h-bonds (LINCS)
Timestep: 2 fs
```

#### Production
```
Duration: 200-500 ns
Temperature: 300 K
Pressure: 1 bar
Barostat: Parrinello-Rahman
Thermostat: V-rescale
Timestep: 2 fs
Output: Every 10 ps (20,000-50,000 frames)
Replica: 3회 (다른 random seed)
```

### 2.3 비결합 상호작용

```
Cutoff scheme: Verlet
vdW cutoff: 1.2 nm
vdW modifier: Force-switch (1.0-1.2 nm)
Coulomb cutoff: 1.2 nm
Coulomb type: PME
PME order: 4
Fourier spacing: 0.12 nm
```

---

## 3. 분석 프로토콜

### 3.1 기본 품질 확인

#### 에너지 안정성
```bash
# Total energy, potential, kinetic
# Temperature, pressure
# 평형 도달 확인
```

**기준:**
- Temperature: 300 ± 5 K
- Pressure: 1 ± 50 bar (fluctuation 정상)
- Energy drift: < 0.1% per ns

#### 구조 안정성
```bash
# RMSD (전체 분자)
# RMSD (Glucose만)
# RMSD (PEG backbone)
```

**기준:**
- RMSD 평형 도달 (plateau)
- Glucose RMSD < 0.3 nm (구조 유지)

---

### 3.2 핵심 분석 항목

#### A. End-to-End 거리

**정의:**
- 분기점 (또는 캡 C 원자) ↔ Glucose C1 (또는 anomeric O)

**분석:**
```python
# 시간에 따른 거리
# 히스토그램 (bin size: 0.1 nm)
# 평균, 표준편차, 최대/최소
```

**예상 분포:**
- 주 피크: 2-3 nm (coil)
- 꼬리: 6-9 nm (extended)
- Bimodal 또는 broad distribution

**Python 예시:**
```python
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

u = mda.Universe('topology.pdb', 'trajectory.xtc')

# 원자 선택
cap_atom = u.select_atoms('name C1 and resid 1')  # 분기점
glc_atom = u.select_atoms('name C1 and resname GLC')  # Glucose

distances = []
for ts in u.trajectory:
    dist = np.linalg.norm(cap_atom.positions - glc_atom.positions)
    distances.append(dist / 10)  # Å to nm

# 히스토그램
plt.hist(distances, bins=50, density=True)
plt.xlabel('End-to-End Distance (nm)')
plt.ylabel('Probability Density')
plt.savefig('end_to_end_distribution.png')

# 통계
print(f"Mean: {np.mean(distances):.2f} nm")
print(f"Std: {np.std(distances):.2f} nm")
print(f"P(R > 6 nm): {np.sum(np.array(distances) > 6) / len(distances):.3f}")
print(f"P(R > 8 nm): {np.sum(np.array(distances) > 8) / len(distances):.3f}")
```

---

#### B. Radius of Gyration (Rg)

**정의:**
- 전체 분자의 질량 중심으로부터 평균 거리

**분석:**
```python
# 시간에 따른 Rg
# 히스토그램
# 평균 ± 표준편차
```

**예상 값:**
- EO24: 1.5-2.5 nm (문헌 기준)

**Python 예시:**
```python
from MDAnalysis.analysis import rg

# 전체 분자 선택
molecule = u.select_atoms('all')

rg_values = []
for ts in u.trajectory:
    rg_val = molecule.radius_of_gyration() / 10  # Å to nm
    rg_values.append(rg_val)

plt.plot(rg_values)
plt.xlabel('Frame')
plt.ylabel('Rg (nm)')
plt.savefig('rg_timeseries.png')

print(f"Rg mean: {np.mean(rg_values):.2f} nm")
print(f"Rg std: {np.std(rg_values):.2f} nm")
```

---

#### C. Tail 확률

**정의:**
- 긴 형태(extended conformation) 출현 빈도

**계산:**
```python
# P(R > 6 nm): 중간 정도 펼쳐짐
# P(R > 8 nm): 거의 완전히 펼쳐짐
```

**해석:**
- P(R > 6 nm) > 0.1: 충분한 펼쳐짐
- P(R > 8 nm) > 0.01: tail 존재 확인

---

#### D. Autocorrelation Function

**목적:**
- 형태 변화 속도 확인
- 평형 도달 검증

**계산:**
```python
from scipy import signal

# End-to-end 거리의 autocorrelation
acf = signal.correlate(distances, distances, mode='full')
acf = acf[len(acf)//2:]
acf /= acf[0]

# Decorrelation time
tau = np.where(acf < 1/np.e)[0][0] * 0.01  # ns (10 ps frame)

print(f"Decorrelation time: {tau:.1f} ns")
```

**기준:**
- tau < 50 ns: 충분한 샘플링 (200 ns에서)
- tau > 100 ns: 시뮬레이션 시간 증가 필요

---

#### E. Dihedral 분포

**중요 Dihedral:**
1. PEG backbone: C-C-O-C
2. Triazole 연결부
3. Amide 연결부
4. O-glycosidic bond

**분석:**
```python
# 각 dihedral의 분포
# Ramachandran-style plot
# 주요 conformer 확인
```

**예상:**
- PEG: gauche (±60°) 선호
- Trans (180°) 일부 존재
- Cis (0°) 거의 없음

---

### 3.3 시각화

#### 필수 그림

1. **End-to-End 거리 히스토그램**
   - 3 replica 겹쳐서 표시
   - 문헌 데이터와 비교

2. **Rg 시계열**
   - 평형 도달 확인
   - 3 replica 비교

3. **대표 구조 스냅샷**
   - Coiled state (R ~ 2-3 nm)
   - Extended state (R ~ 7-9 nm)
   - VMD/PyMOL 렌더링

4. **Dihedral 분포**
   - PEG backbone
   - 연결부

---

## 4. 판정 기준

### ✅ 합격 (CGenFF 사용 가능)

- [ ] End-to-end 평균: 2-4 nm
- [ ] P(R > 6 nm): 0.05-0.20
- [ ] P(R > 8 nm): 0.001-0.05
- [ ] Rg: 1.5-2.5 nm
- [ ] Autocorrelation time < 50 ns
- [ ] 3 replica 일관성 (평균 ± 20%)
- [ ] 에너지 안정
- [ ] Glucose 구조 유지 (RMSD < 0.3 nm)

### ⚠️ 재검토 필요

**과도한 뭉침:**
- End-to-end < 2 nm (대부분)
- P(R > 6 nm) < 0.01
- Rg < 1.2 nm

**비현실적 펼쳐짐:**
- End-to-end > 5 nm (대부분)
- P(R > 8 nm) > 0.2
- Rg > 3.5 nm

**동역학 문제:**
- Autocorrelation time > 100 ns
- 평형 미도달
- Replica 간 큰 차이 (> 50%)

---

## 5. 문제 해결

### 시나리오 1: 과도한 뭉침

**진단:**
```bash
# 1. CGenFF penalty 재확인
grep "PENALTY" ligand.str | sort -k2 -n -r | head -20

# 2. 비결합 파라미터 확인
# LJ epsilon, sigma 값
# 전하 분포

# 3. Dihedral 분포 확인
# 특정 dihedral이 한 곳에 갇혀있는지
```

**해결:**
1. High penalty dihedral 재파라미터화
2. Enhanced sampling (replica exchange)
3. 시뮬레이션 시간 증가 (500 ns - 1 μs)

---

### 시나리오 2: Replica 간 불일치

**원인:**
- 샘플링 부족
- 에너지 장벽 높음

**해결:**
1. 시뮬레이션 시간 증가
2. Replica 수 증가 (5-10개)
3. 다른 초기 구조 사용

---

### 시나리오 3: Glucose 구조 변형

**원인:**
- CGenFF Glucose 파라미터 문제
- 연결부 strain

**해결:**
1. CHARMM36 당 파라미터 사용 고려
2. O-glycosidic bond 파라미터 재검토
3. Restraint 사용 (임시)

---

## 6. 보고서 작성

### 필수 포함 항목

1. **시스템 정보**
   - 분자 구조 (SMILES, 그림)
   - Force field (CGenFF version)
   - Penalty score 요약

2. **시뮬레이션 조건**
   - Box size, 이온 농도
   - 온도, 압력
   - 시뮬레이션 시간, replica 수

3. **결과**
   - End-to-end 분포 (그림 + 통계)
   - Rg 분포
   - Tail 확률
   - Autocorrelation time

4. **대표 구조**
   - Coiled, extended 스냅샷
   - VMD 렌더링

5. **결론**
   - CGenFF 사용 가능 여부
   - 3-arm 진행 권장 사항
   - 추가 최적화 필요 항목

---

## 7. 다음 단계

### CGenFF 합격 시

1. 3-arm 시스템 구축
2. 동일 조건 MD (더 긴 시간)
3. 팔-팔 상호작용 분석

### 재검토 필요 시

1. Dihedral 재파라미터화
2. 1-arm 재검증
3. 대안 force field 고려 (GLYCAM 등)

---

**작성자:** Cascade AI  
**최종 수정:** 2026-01-18
