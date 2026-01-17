# Glycocalyx 모델 개선 가이드
# 현재 설정 분석 및 개선 방안

## 📊 현재 설정 분석

### ✅ 잘 설정된 부분

**Glycan 구조**:
```
aDMan(1→3)[aDMan(1→6)]bDMan(1→4)bDGlcNAc(1→4)bDGlcNAc(1→)PROA-45
```

- **Type**: High-Mannose N-glycan
- **위치**: Asn45 (GLUT1의 extracellular loop)
- **구성**: 2 GlcNAc + 3 Mannose (총 5개 당)
- **크기**: 약 2-3 nm

**장점**:
- ✅ N-glycan의 기본 구조 포함
- ✅ Asn45에 정확히 부착
- ✅ 입체 장벽 효과 일부 제공

---

## ⚠️ 주요 개선 필요 사항

### 1. 시알산(Sialic Acid) 부재 ⭐⭐⭐ (최우선)

**현재 문제**:
- 현재 glycan은 **중성 당**만 포함 (Mannose, GlcNAc)
- **음전하 없음** → Tripod와의 정전기적 반발 효과 없음
- RBC glycocalyx의 핵심 특성 미반영

**RBC Glycocalyx의 핵심**:
- **시알산 (NeuAc, Neu5Ac)**: 음전하 (-1) 제공
- RBC 표면 전하 밀도: -0.02 to -0.05 C/m²
- Tripod의 L-glucose (중성/약양성)와 정전기적 반발

**해결 방안**:
→ **Complex Type N-glycan**으로 변경하고 **시알산 추가**

---

### 2. Glycan 크기 부족

**현재**: 5개 당 (2-3 nm)
**필요**: RBC glycocalyx는 5-10 nm 두께

**문제점**:
- Tripod (직경 약 3-4 nm)가 쉽게 통과 가능
- 입체 장벽 효과 불충분

**해결 방안**:
→ 더 긴 Complex Type N-glycan 사용

---

### 3. Glycan 개수 부족

**현재**: Asn45 1개 부위만
**GLUT1의 N-glycosylation sites**: 14개 ASN 잔기

**주요 extracellular sites**:
- Asn45 ✅ (현재 설정됨)
- Asn88
- Asn100
- Asn288
- Asn317

**문제점**:
- 1개 glycan으로는 glycocalyx "숲" 효과 부족
- RBC는 표면이 glycan으로 빽빽하게 덮여 있음

**해결 방안**:
→ 최소 3-5개 부위에 glycan 추가

---

## 🛠️ 구체적 개선 방안

### 개선안 A: Complex Type N-glycan with Sialic Acid (추천)

#### CHARMM-GUI 설정

**Step 1: Glycan 구조 변경**

현재 (High-Mannose):
```
Man3GlcNAc2
```

변경 후 (Complex Type with Sialic Acid):
```
Neu5Ac(α2-3)Gal(β1-4)GlcNAc(β1-2)Man(α1-3)
                                        \
                                     GlcNAc(β1-4)GlcNAc(β1-)Asn
                                        /
Neu5Ac(α2-3)Gal(β1-4)GlcNAc(β1-2)Man(α1-6)
```

**구성**:
- Core: GlcNAc-GlcNAc-Man (동일)
- Branches: 2개 (bi-antennary)
- 각 branch: Man → GlcNAc → Gal → **Neu5Ac (시알산)**
- 총 당 개수: 9개 (크기 약 4-5 nm)
- **음전하**: -2 (시알산 2개)

#### CHARMM-GUI 입력 방법

**Option 1: Glycan Reader - Manual Input**

```
Glycan Sequence Builder:

Branch 1:
β-D-GlcNAc(1→4)β-D-GlcNAc(1→)Asn45
    ↓
α-D-Man(1→3)
    ↓
β-D-GlcNAc(1→2)
    ↓
β-D-Gal(1→4)
    ↓
α-Neu5Ac(2→3)  ← 시알산 추가!

Branch 2:
α-D-Man(1→6) [from core Man]
    ↓
β-D-GlcNAc(1→2)
    ↓
β-D-Gal(1→4)
    ↓
α-Neu5Ac(2→3)  ← 시알산 추가!
```

**Option 2: Glycan Reader - Preset Selection**

CHARMM-GUI에서 제공하는 preset 사용:
1. **"Complex Type N-glycan"** 선택
2. **"Bi-antennary"** 선택
3. **"Add terminal sialic acid (α2-3 linkage)"** 체크
4. Attachment site: **Asn45**

#### GRS (Glycan Representation Sequence) 코드

CHARMM-GUI Upload GRS 기능 사용 시:

```
BGLCNA(1-4)BGLCNA(1-N)PROA-45
AMAN(1-3)[AMAN(1-6)]BMAN(1-4)BGLCNA
BGLCNA(1-2)AMAN(1-3)BMAN
BGAL(1-4)BGLCNA
ANE5AC(2-3)BGAL
BGLCNA(1-2)AMAN(1-6)BMAN
BGAL(1-4)BGLCNA
ANE5AC(2-3)BGAL
```

---

### 개선안 B: 여러 Glycosylation Site 활용

#### GLUT1의 주요 N-glycosylation sites

| Site | 위치 | 접근성 | 우선순위 |
|------|------|--------|----------|
| **Asn45** | Extracellular loop 1 | 높음 | ⭐⭐⭐ (현재 설정) |
| **Asn88** | Extracellular loop 2 | 높음 | ⭐⭐⭐ |
| **Asn100** | Extracellular loop 2 | 중간 | ⭐⭐ |
| **Asn288** | Extracellular loop 7 | 높음 | ⭐⭐⭐ |
| **Asn317** | Extracellular loop 8 | 중간 | ⭐⭐ |

#### 추천 설정

**Minimal Setup (3 sites)**:
- Asn45: Complex Type + Sialic Acid
- Asn88: Complex Type + Sialic Acid
- Asn288: Complex Type + Sialic Acid

**Optimal Setup (5 sites)**:
- Asn45, 88, 100, 288, 317: 모두 Complex Type + Sialic Acid

#### CHARMM-GUI 설정 방법

1. **Glycan Reader** 재실행
2. **Multiple Glycosylation Sites** 선택
3. 각 site에 동일한 glycan 구조 적용:
   ```
   Site 1: Asn45 → Complex bi-antennary + Neu5Ac
   Site 2: Asn88 → Complex bi-antennary + Neu5Ac
   Site 3: Asn288 → Complex bi-antennary + Neu5Ac
   ```

---

### 개선안 C: Glycolipid 추가 (보완)

Glycoprotein만으로 부족할 경우, **Glycolipid (GM1, GM3)**를 막에 추가:

#### CHARMM-GUI Membrane Builder 설정

**Lipid Composition (Outer Leaflet)**:
```
- PSM: 22%
- POPC: 18%
- GM1: 8%    ← Ganglioside (시알산 1개)
- GM3: 5%    ← Ganglioside (시알산 1개)
- CHL1: 47%
```

**GM1 구조**:
```
Neu5Ac(α2-3)Gal(β1-3)GalNAc(β1-4)[Gal(β1-3)]Glc(β1-1)Ceramide
```
- 시알산 1개 포함
- 음전하 -1
- 막에 직접 삽입

**장점**:
- Glycoprotein + Glycolipid = 더 조밀한 glycocalyx
- 막 표면 전체에 음전하 분포
- RBC 환경에 더 가까움

---

## 📋 단계별 실행 계획

### Phase 1: 시알산 추가 (최우선) ⭐⭐⭐

**목표**: 현재 Asn45 glycan에 시알산 추가

**방법**:
1. CHARMM-GUI Glycan Reader 재실행
2. 기존 PDB 업로드
3. Asn45의 glycan을 **Complex Type + Sialic Acid**로 변경
4. 새로운 시스템 다운로드

**예상 효과**:
- 음전하 -2 추가
- Tripod와의 정전기적 반발 발생
- Glycocalyx 장벽 효과 증가

**검증 방법**:
```python
# Electrostatic potential 계산
from pdbfixer import PDBFixer
import openmm as mm

# 시알산의 전하 확인
# COO- 그룹이 -1 전하를 띠는지 확인
```

---

### Phase 2: Glycan 개수 증가 ⭐⭐

**목표**: 3-5개 부위에 glycan 추가

**방법**:
1. CHARMM-GUI Glycan Reader
2. Multiple sites 선택:
   - Asn45, Asn88, Asn288 (최소)
   - + Asn100, Asn317 (최적)
3. 모든 site에 동일한 Complex Type glycan 적용

**예상 효과**:
- Glycocalyx 밀도 증가
- 입체 장벽 효과 강화
- Tripod 접근 차단 확률 증가

---

### Phase 3: Glycolipid 추가 (선택) ⭐

**목표**: 막에 GM1/GM3 ganglioside 추가

**방법**:
1. CHARMM-GUI Membrane Builder
2. Lipid composition 수정
3. GM1 8%, GM3 5% 추가

**예상 효과**:
- 막 표면 전체에 음전하 분포
- 더 현실적인 RBC 환경
- Glycocalyx 밀도 최대화

---

## 🧪 시뮬레이션 및 검증

### 검증 항목

#### 1. Electrostatic Potential Map

```python
# scripts/analyze_electrostatic.py
import mdtraj as md
import numpy as np
from gridData import Grid

def calculate_surface_charge(pdb_file):
    """
    Glycocalyx 표면의 전하 분포 계산
    """
    traj = md.load(pdb_file)
    
    # 시알산 잔기 선택
    sialic_acids = traj.top.select('resname ANE5AC or resname BNEN')
    
    # 전하 계산
    charges = []
    for atom in sialic_acids:
        if 'COO' in atom.name:  # Carboxyl group
            charges.append(-1.0)
    
    total_charge = sum(charges)
    print(f"Total negative charge from sialic acids: {total_charge}")
    
    return total_charge

# 실행
charge_rbc = calculate_surface_charge('rbc_glycocalyx.pdb')
charge_tumor = calculate_surface_charge('tumor_control.pdb')

print(f"\nRBC surface charge: {charge_rbc}")
print(f"Tumor surface charge: {charge_tumor}")
print(f"Charge difference: {charge_rbc - charge_tumor}")
```

**예상 결과**:
- RBC (with sialic acid): -6 to -10 (3-5 glycans × 2 sialic acids)
- Tumor (no glycocalyx): 0 to -2
- **차이가 클수록 선택성 증가**

#### 2. Glycocalyx Thickness

```python
# scripts/measure_glycocalyx_thickness.py
import mdtraj as md
import numpy as np

def measure_glycocalyx_thickness(traj):
    """
    Glycocalyx 층의 두께 측정
    """
    # 막 표면 (phosphate 그룹)
    membrane_surface = traj.atom_slice(
        traj.top.select('name P')
    )
    membrane_z = membrane_surface.xyz[:, :, 2].mean()
    
    # Glycan 끝 (시알산)
    glycan_tips = traj.atom_slice(
        traj.top.select('resname ANE5AC')
    )
    glycan_z_max = glycan_tips.xyz[:, :, 2].max()
    
    thickness = (glycan_z_max - membrane_z) * 10  # nm to Å
    
    return thickness

# 실행
thickness = measure_glycocalyx_thickness(traj)
print(f"Glycocalyx thickness: {thickness:.2f} Å")
```

**목표**:
- High-Mannose (현재): 20-30 Å
- Complex + Sialic Acid: **40-50 Å** (목표)
- RBC 실제: 50-100 Å

#### 3. Penetration Test

```python
# scripts/test_tripod_penetration.py
import mdtraj as md
import numpy as np

def test_penetration(traj, ligand_name='LIG'):
    """
    Tripod가 glycocalyx를 침투했는지 확인
    """
    # Tripod 위치
    ligand = traj.atom_slice(traj.top.select(f'resname {ligand_name}'))
    ligand_z = ligand.xyz[:, :, 2].mean(axis=1)
    
    # Glycocalyx 경계
    glycan = traj.atom_slice(traj.top.select('resname BGLCNA or resname ANE5AC'))
    glycan_z_min = glycan.xyz[:, :, 2].min(axis=1)
    
    # 침투 여부
    penetration = ligand_z < glycan_z_min
    penetration_ratio = np.sum(penetration) / len(traj)
    
    return penetration_ratio

# 실행
pen_rbc = test_penetration(traj_rbc)
pen_tumor = test_penetration(traj_tumor)

print(f"RBC penetration: {pen_rbc:.2%}")
print(f"Tumor penetration: {pen_tumor:.2%}")
```

**예상 결과**:
- RBC (with glycocalyx): < 10% (차단 성공)
- Tumor (no glycocalyx): > 80% (접근 허용)

---

## 📊 예상 결과 비교

### Before (현재 설정)

| 항목 | 값 | 평가 |
|------|-----|------|
| Glycan type | High-Mannose | ⚠️ 중성 |
| Sialic acid | 0개 | ❌ 없음 |
| Surface charge | 0 | ❌ 중성 |
| Thickness | 20-30 Å | ⚠️ 얇음 |
| Glycan sites | 1개 (Asn45) | ⚠️ 부족 |
| Barrier effect | 약함 | ⚠️ 불충분 |

### After (개선 후)

| 항목 | 값 | 평가 |
|------|-----|------|
| Glycan type | Complex bi-antennary | ✅ 현실적 |
| Sialic acid | 6-10개 | ✅ 충분 |
| Surface charge | -6 to -10 | ✅ 음전하 |
| Thickness | 40-50 Å | ✅ 적절 |
| Glycan sites | 3-5개 | ✅ 충분 |
| Barrier effect | 강함 | ✅ 효과적 |

---

## 🎯 최종 추천 설정

### Recommended Setup (균형)

**Glycoprotein**:
- **3개 sites**: Asn45, Asn88, Asn288
- **Glycan type**: Complex bi-antennary
- **Terminal**: α-Neu5Ac(2-3) sialic acid
- **Total charge**: -6 (3 sites × 2 sialic acids)

**Membrane** (선택적):
- GM1: 5-8%
- GM3: 3-5%
- 나머지: 기존 조성 유지

**예상 효과**:
- ✅ 충분한 음전하 (-6 to -10)
- ✅ 적절한 glycocalyx 두께 (40-50 Å)
- ✅ 입체 장벽 효과
- ✅ Tripod 선택성 검증 가능

---

## 💡 추가 팁

### Tip 1: 시알산 Linkage 선택

**α(2-3) vs α(2-6)**:
- **α(2-3)**: RBC에서 더 흔함 (추천)
- **α(2-6)**: 암세포에서 증가

→ RBC 모델에는 **α(2-3)** 사용

### Tip 2: Glycan 다양성

실제 RBC는 다양한 glycan 구조를 가짐:
- 일부는 High-Mannose
- 일부는 Complex Type
- 일부는 Hybrid Type

→ 더 현실적으로 만들려면 **혼합 사용**

### Tip 3: Force Field 확인

CHARMM36 carbohydrate force field가 시알산을 지원하는지 확인:
```bash
# toppar 폴더 확인
grep -r "NEU5AC\|ANE5AC" toppar/
```

---

## ✅ Action Items

### 즉시 실행 (Priority 1)

1. [ ] CHARMM-GUI Glycan Reader 재실행
2. [ ] Asn45 glycan을 Complex Type + Sialic Acid로 변경
3. [ ] 새 시스템 다운로드 및 검증

### 단기 실행 (Priority 2)

4. [ ] Asn88, Asn288에 glycan 추가
5. [ ] Electrostatic potential map 생성
6. [ ] Glycocalyx thickness 측정

### 장기 실행 (Priority 3)

7. [ ] Glycolipid (GM1/GM3) 추가
8. [ ] 다양한 glycan 구조 테스트
9. [ ] 최적 조합 도출

---

**질문이나 추가 설명이 필요하면 알려주세요!**
