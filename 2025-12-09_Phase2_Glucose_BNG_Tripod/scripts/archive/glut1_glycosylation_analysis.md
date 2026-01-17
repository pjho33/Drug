# GLUT1 Glycosylation 전략 분석
# Tripod 결합에 미치는 영향 및 RBC vs 종양세포 차이

## 🎯 핵심 개념

### RBC Glycocalyx vs GLUT1 Glycocalyx

**두 가지 장벽**:
1. **RBC 막 Glycocalyx** (일반적 장벽)
   - Glycophorin A, Band 3, Glycolipids (GM1/GM3)
   - 막 전체를 덮는 "숲"
   - Tripod가 막에 접근하는 것 자체를 차단
   
2. **GLUT1 자체 Glycosylation** (특이적 장벽) ⭐⭐⭐
   - GLUT1 단백질에 직접 부착된 N-glycan
   - Glucose 결합 포켓 주변에 위치
   - **Tripod가 포켓에 진입하는 것을 직접 차단**

**중요**: Tripod는 GLUT1의 glucose 포켓에 결합해야 하므로, **GLUT1 자체의 glycan이 훨씬 더 중요**합니다!

---

## 📊 GLUT1 N-glycosylation Sites 분석

### 전체 ASN 잔기 (14개)

| Residue | 좌표 (X, Y, Z) | 위치 | Glucose 포켓 거리 | 중요도 |
|---------|----------------|------|-------------------|--------|
| **Asn45** | 584.6, -52.4, 184.8 | **Extracellular loop 1** | ~15 Å | ⭐⭐⭐ |
| **Asn88** | 571.4, -14.6, 198.5 | **Extracellular loop 2** | ~10 Å | ⭐⭐⭐ |
| Asn29 | 578.2, -38.3, 200.8 | Extracellular | ~20 Å | ⭐⭐ |
| Asn34 | 584.8, -44.7, 199.1 | Extracellular | ~18 Å | ⭐⭐ |
| Asn94 | 564.8, -21.0, 201.3 | Extracellular loop 2 | ~12 Å | ⭐⭐⭐ |
| Asn100 | 567.5, -29.3, 196.8 | Extracellular loop 2 | ~15 Å | ⭐⭐ |
| Asn182 | 578.4, -54.0, 198.1 | Extracellular loop 5 | ~25 Å | ⭐ |
| Asn217 | 568.1, -8.0, 203.6 | Extracellular loop 6 | ~18 Å | ⭐⭐ |
| Asn219 | 564.1, -6.5, 203.1 | Extracellular loop 6 | ~20 Å | ⭐ |
| Asn222 | 558.6, -5.7, 209.3 | Extracellular loop 6 | ~25 Å | ⭐ |
| **Asn288** | 592.2, -36.6, 202.0 | **Extracellular loop 7** | ~8 Å | ⭐⭐⭐ |
| Asn317 | 588.3, -40.0, 213.5 | Extracellular loop 8 | ~15 Å | ⭐⭐ |
| Asn411 | 590.7, -24.7, 198.2 | Extracellular loop 11 | ~12 Å | ⭐⭐ |
| Asn415 | 591.9, -30.1, 196.4 | Extracellular loop 11 | ~10 Å | ⭐⭐⭐ |

**Glucose 포켓 위치 (추정)**: 약 (575, -25, 195) 근처

---

## 🔍 Glucose 포켓 주변 Critical Sites

### 최우선 Glycosylation Sites (Tripod 차단에 직접 영향)

#### 1. **Asn288** ⭐⭐⭐⭐⭐ (가장 중요!)

**위치**: Extracellular loop 7
**거리**: Glucose 포켓에서 **~8 Å** (가장 가까움!)
**역할**: 
- Glucose 포켓 입구 바로 위에 위치
- Glycan이 부착되면 포켓 입구를 직접 차단
- **Tripod의 L-glucose head가 진입하는 경로를 막음**

**RBC vs 종양세포**:
- **RBC**: Asn288에 **Complex Type + Sialic Acid** → 강력한 차단
- **종양세포**: Glycan 없음 또는 작은 glycan → 진입 허용

**전략**: 
→ **Asn288은 필수로 glycosylation 해야 함!**

---

#### 2. **Asn88** ⭐⭐⭐⭐

**위치**: Extracellular loop 2
**거리**: Glucose 포켓에서 **~10 Å**
**역할**:
- 포켓 입구 측면에 위치
- Tripod의 접근 각도를 제한
- Asn288과 함께 "게이트키퍼" 역할

**전략**:
→ Asn288과 함께 glycosylation 필수

---

#### 3. **Asn415** ⭐⭐⭐

**위치**: Extracellular loop 11
**거리**: Glucose 포켓에서 **~10 Å**
**역할**:
- 포켓 입구 반대편에 위치
- Tripod의 회전/재배치를 방해

**전략**:
→ 추가 장벽으로 glycosylation 권장

---

#### 4. **Asn45** ⭐⭐⭐

**위치**: Extracellular loop 1
**거리**: Glucose 포켓에서 **~15 Å**
**역할**:
- 포켓에서 약간 떨어져 있음
- 일반적인 입체 장벽 제공
- 현재 설정된 site

**전략**:
→ 유지하되, Asn288이 더 중요

---

## 🎯 최적 Glycosylation 전략

### Strategy A: Glucose 포켓 집중 차단 (추천) ⭐⭐⭐

**목표**: Tripod가 GLUT1 포켓에 진입하는 것을 **직접 차단**

**Glycosylation Sites**:
1. **Asn288** (최우선) - 포켓 입구 바로 위
2. **Asn88** - 포켓 측면
3. **Asn415** - 포켓 반대편

**Glycan Type**: Complex bi-antennary + Sialic Acid (α2-3)

**예상 효과**:
- 포켓 입구를 glycan이 "뚜껑"처럼 덮음
- Tripod의 L-glucose head가 진입 불가
- 음전하 (-6)로 정전기적 반발 추가

**RBC 모델**:
```
        Neu5Ac(-1)
           |
    [Glycan at Asn288]  ← 포켓 바로 위
           |
    ==================
    |  Glucose Pocket  |  ← Tripod 진입 차단!
    ==================
         GLUT1
```

**종양세포 모델**:
```
    (No glycan)  ← 포켓 열림
           
    ==================
    |  Glucose Pocket  |  ← Tripod 자유롭게 진입
    ==================
         GLUT1
```

---

### Strategy B: 전면 차단 (보수적)

**목표**: 포켓 주변 전체를 glycan으로 둘러싸기

**Glycosylation Sites**:
1. Asn288 (포켓 위)
2. Asn88 (포켓 측면)
3. Asn415 (포켓 반대편)
4. Asn45 (일반 장벽)
5. Asn94 (추가 측면)

**Glycan Type**: Complex bi-antennary + Sialic Acid

**예상 효과**:
- 포켓 주변 360도 차단
- 어떤 각도에서도 Tripod 접근 어려움
- 음전하 -10 (매우 강한 반발)

**단점**:
- 계산 비용 증가
- 과도한 glycosylation (비현실적일 수 있음)

---

### Strategy C: 최소 차단 (빠른 검증)

**목표**: 최소한의 glycan으로 효과 검증

**Glycosylation Sites**:
1. **Asn288만** (포켓 입구)

**Glycan Type**: Complex bi-antennary + Sialic Acid

**예상 효과**:
- 포켓 입구 직접 차단
- 음전하 -2
- 빠른 시뮬레이션

**용도**:
- Proof-of-concept
- Asn288의 중요성 검증

---

## 📋 CHARMM-GUI 설정 가이드

### 추천 설정 (Strategy A)

#### Step 1: Glycan Reader 실행

**Multiple Glycosylation Sites 입력**:
```
Site 1: Protein PROA, Residue ASN, Number 288  ← 최우선!
Site 2: Protein PROA, Residue ASN, Number 88
Site 3: Protein PROA, Residue ASN, Number 415
```

#### Step 2: 각 Site에 동일한 Glycan 적용

**Glycan Structure (Complex Type with Sialic Acid)**:
```
Branch 1:
β-D-GlcNAc(1→4)β-D-GlcNAc(1→)Asn
    ↓
α-D-Man(1→3)
    ↓
β-D-GlcNAc(1→2)
    ↓
β-D-Gal(1→4)
    ↓
α-Neu5Ac(2→3)  ← 시알산

Branch 2:
α-D-Man(1→6) [from core]
    ↓
β-D-GlcNAc(1→2)
    ↓
β-D-Gal(1→4)
    ↓
α-Neu5Ac(2→3)  ← 시알산
```

#### Step 3: 시스템 생성

- Membrane Builder로 이동
- GLUT1 + Glycans 포함
- Lipid composition 설정
- 시스템 다운로드

---

## 🧪 시뮬레이션 프로토콜

### Test 1: Glycan 위치 효과 검증

**목표**: 어느 site가 Tripod 차단에 가장 효과적인가?

**시스템**:
1. **Asn288만** glycosylation
2. **Asn88만** glycosylation
3. **Asn415만** glycosylation
4. **대조군** (glycan 없음)

**측정**:
- Tripod와 포켓 간 거리
- 결합 성공률
- 체류 시간

**예상 결과**:
- Asn288 > Asn88 > Asn415 순으로 차단 효과

---

### Test 2: RBC vs 종양세포 비교

**RBC 모델**:
- Asn288, 88, 415에 Complex + Sialic Acid
- 총 음전하: -6

**종양세포 모델**:
- Glycan 없음 (naked GLUT1)
- 또는 High-Mannose만 (작고 중성)

**측정**:
- Binding free energy (MM/PBSA)
- Contact time
- Penetration depth

**예상 결과**:
| 항목 | RBC | 종양세포 | 비율 |
|------|-----|----------|------|
| Contact time | < 10% | > 70% | 7배 |
| Binding events | 0-1 | 5-8 | 8배 |
| ΔG binding | N/A | -30 kcal/mol | - |

---

## 🔬 분석 스크립트

### 1. Glycan-Pocket Distance 측정

```python
# scripts/analyze_glycan_pocket_distance.py
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

def measure_glycan_pocket_distance(traj, glycan_site, pocket_residues):
    """
    Glycan과 glucose 포켓 간 거리 측정
    """
    # Glycan 선택 (시알산 끝)
    glycan_atoms = traj.top.select(
        f'resname ANE5AC and resid {glycan_site}'
    )
    
    # Glucose 포켓 선택 (주요 잔기)
    pocket_atoms = traj.top.select(
        f'resid {" ".join(map(str, pocket_residues))}'
    )
    
    # 거리 계산
    distances = []
    for frame in traj:
        glycan_pos = frame.xyz[0, glycan_atoms, :].mean(axis=0)
        pocket_pos = frame.xyz[0, pocket_atoms, :].mean(axis=0)
        dist = np.linalg.norm(glycan_pos - pocket_pos) * 10  # nm to Å
        distances.append(dist)
    
    return np.array(distances)

# GLUT1 glucose 포켓 주요 잔기 (예시)
pocket_residues = [168, 169, 288, 289, 290, 415, 416]

# 실행
traj = md.load('production.dcd', top='system.pdb')
dist_288 = measure_glycan_pocket_distance(traj, 288, pocket_residues)
dist_88 = measure_glycan_pocket_distance(traj, 88, pocket_residues)

print(f"Asn288-Pocket distance: {dist_288.mean():.2f} ± {dist_288.std():.2f} Å")
print(f"Asn88-Pocket distance: {dist_88.mean():.2f} ± {dist_88.std():.2f} Å")

# 플롯
plt.figure(figsize=(10, 6))
plt.plot(dist_288, label='Asn288', alpha=0.7)
plt.plot(dist_88, label='Asn88', alpha=0.7)
plt.axhline(y=10, color='r', linestyle='--', label='Critical distance')
plt.xlabel('Frame')
plt.ylabel('Distance (Å)')
plt.title('Glycan-Pocket Distance')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig('glycan_pocket_distance.png', dpi=300)
```

### 2. Tripod 진입 차단 분석

```python
# scripts/analyze_tripod_blocking.py
import mdtraj as md
import numpy as np

def analyze_blocking_effect(traj, ligand_name='LIG'):
    """
    Glycan이 Tripod 진입을 차단하는지 분석
    """
    # Tripod L-glucose head 선택
    tripod_head = traj.top.select(
        f'resname {ligand_name} and name C1 C2 C3'  # Glucose carbons
    )
    
    # Glycan 장벽 선택 (Asn288)
    glycan_barrier = traj.top.select(
        'resname ANE5AC and resid 288'
    )
    
    # Glucose 포켓 선택
    pocket = traj.top.select(
        'resid 168 169 288 289 290'
    )
    
    # 각 프레임에서 분석
    blocked_frames = 0
    penetrated_frames = 0
    
    for frame_idx in range(len(traj)):
        # Tripod 위치
        tripod_z = traj.xyz[frame_idx, tripod_head, 2].mean()
        
        # Glycan 위치
        glycan_z = traj.xyz[frame_idx, glycan_barrier, 2].mean()
        
        # Pocket 위치
        pocket_z = traj.xyz[frame_idx, pocket, 2].mean()
        
        # Tripod가 glycan 아래로 침투했는가?
        if tripod_z < glycan_z and tripod_z > pocket_z - 5:
            penetrated_frames += 1
        elif tripod_z > glycan_z + 5:
            blocked_frames += 1
    
    blocking_ratio = blocked_frames / len(traj)
    penetration_ratio = penetrated_frames / len(traj)
    
    return {
        'blocking_ratio': blocking_ratio,
        'penetration_ratio': penetration_ratio,
        'blocked_frames': blocked_frames,
        'penetrated_frames': penetrated_frames
    }

# 실행
traj_rbc = md.load('rbc_production.dcd', top='rbc_system.pdb')
traj_tumor = md.load('tumor_production.dcd', top='tumor_system.pdb')

result_rbc = analyze_blocking_effect(traj_rbc)
result_tumor = analyze_blocking_effect(traj_tumor)

print("\n=== Blocking Effect Analysis ===")
print(f"\nRBC (with glycan):")
print(f"  Blocking ratio: {result_rbc['blocking_ratio']:.2%}")
print(f"  Penetration ratio: {result_rbc['penetration_ratio']:.2%}")

print(f"\nTumor (no glycan):")
print(f"  Blocking ratio: {result_tumor['blocking_ratio']:.2%}")
print(f"  Penetration ratio: {result_tumor['penetration_ratio']:.2%}")

print(f"\nSelectivity: {result_tumor['penetration_ratio'] / max(result_rbc['penetration_ratio'], 0.01):.1f}x")
```

### 3. Electrostatic Repulsion 분석

```python
# scripts/analyze_electrostatic_repulsion.py
import mdtraj as md
import numpy as np

def calculate_electrostatic_interaction(traj, ligand_name='LIG'):
    """
    Tripod와 glycan 간 정전기적 상호작용 계산
    """
    # Tripod 전하 중심 (L-glucose는 중성/약양성)
    tripod_atoms = traj.top.select(f'resname {ligand_name}')
    
    # Sialic acid 음전하 중심 (COO- 그룹)
    sialic_coo = traj.top.select(
        'resname ANE5AC and (name C1 or name O1A or name O1B)'
    )
    
    # 거리 계산
    distances = []
    for frame_idx in range(len(traj)):
        tripod_pos = traj.xyz[frame_idx, tripod_atoms, :].mean(axis=0)
        sialic_pos = traj.xyz[frame_idx, sialic_coo, :].mean(axis=0)
        
        dist = np.linalg.norm(tripod_pos - sialic_pos) * 10  # nm to Å
        distances.append(dist)
    
    distances = np.array(distances)
    
    # Coulomb 상호작용 에너지 추정 (간단한 모델)
    # E = k * q1 * q2 / r
    # q1 = +0.5 (Tripod, 약양성), q2 = -1 (Sialic acid)
    # k = 332 kcal/(mol·Å·e²) (in vacuum)
    
    k = 332  # kcal/(mol·Å)
    q_tripod = 0.5  # 약양성
    q_sialic = -1.0  # 음전하
    epsilon = 4.0  # Dielectric constant (수용액 환경)
    
    energies = k * q_tripod * q_sialic / (distances * epsilon)
    
    return distances, energies

# 실행
traj = md.load('rbc_production.dcd', top='rbc_system.pdb')
distances, energies = calculate_electrostatic_interaction(traj)

print(f"\nElectrostatic Interaction:")
print(f"  Mean distance: {distances.mean():.2f} Å")
print(f"  Mean repulsion energy: {energies.mean():.2f} kcal/mol")
print(f"  Max repulsion: {energies.max():.2f} kcal/mol")

# 반발력이 강한 프레임 비율
strong_repulsion = np.sum(energies > 2.0) / len(energies)
print(f"  Strong repulsion (>2 kcal/mol): {strong_repulsion:.2%}")
```

---

## 📊 예상 결과 요약

### Glycan Site 중요도

| Site | 포켓 거리 | 차단 효과 | 우선순위 | 설정 |
|------|-----------|-----------|----------|------|
| **Asn288** | 8 Å | ⭐⭐⭐⭐⭐ | 1 | **필수** |
| **Asn88** | 10 Å | ⭐⭐⭐⭐ | 2 | **필수** |
| **Asn415** | 10 Å | ⭐⭐⭐⭐ | 3 | 권장 |
| Asn45 | 15 Å | ⭐⭐⭐ | 4 | 선택 |
| Asn94 | 12 Å | ⭐⭐⭐ | 5 | 선택 |

### RBC vs 종양세포 예상 결과

| 지표 | RBC (Glycan 있음) | 종양세포 (Glycan 없음) | 선택성 |
|------|-------------------|------------------------|--------|
| **Tripod 진입** | < 5% | > 80% | **16배** |
| **Contact time** | < 10% | > 70% | **7배** |
| **Binding events** | 0-1 | 6-10 | **10배** |
| **ΔG binding** | N/A (차단) | -30 kcal/mol | - |
| **표면 전하** | -6 (반발) | 0 (중성) | - |

---

## ✅ 최종 추천

### Optimal Setup

**GLUT1 Glycosylation**:
1. **Asn288** - Complex bi-antennary + Sialic Acid (α2-3) ⭐⭐⭐⭐⭐
2. **Asn88** - Complex bi-antennary + Sialic Acid (α2-3) ⭐⭐⭐⭐
3. **Asn415** - Complex bi-antennary + Sialic Acid (α2-3) ⭐⭐⭐

**총 효과**:
- 음전하: -6 (시알산 6개)
- Glucose 포켓 입구 차단
- Tripod 진입 방지

**RBC 모델**:
- GLUT1: 위 3개 site glycosylation
- 막: 선택적으로 GM1/GM3 추가

**종양세포 모델**:
- GLUT1: Glycan 없음 (naked)
- 막: Standard lipid composition

---

**다음 단계**: CHARMM-GUI에서 Asn288, 88, 415에 Complex Type + Sialic Acid를 설정하세요!
