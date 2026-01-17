# RBC Glycocalyx 모델링 가이드
# Tripod 약물의 선택적 종양 타겟팅 검증을 위한 시뮬레이션

## 🎯 목표

**가설**: Tripod 약물은 glycocalyx 장벽이 두꺼운 RBC/내피세포에는 접근이 어렵지만, glycocalyx가 불균일한 종양세포의 GLUT1에는 효과적으로 결합할 수 있다.

**검증 방법**: 
1. RBC 모델 (Glycocalyx 있음) - Tripod 접근 차단
2. 종양세포 모델 (Glycocalyx 약함/없음) - Tripod 접근 허용
3. 대조군 (Naked membrane) - 기본 결합력 측정

---

## 📊 Phase 1: RBC Glycocalyx 구성 이해

### 1.1 RBC 막 지질 조성

RBC는 일반 세포와 다른 지질 조성을 가집니다:

**외막 (Outer Leaflet)**
- Sphingomyelin (PSM): 30%
- Phosphatidylcholine (POPC): 25%
- Cholesterol (CHL1): 45% ← **매우 높음!**

**내막 (Inner Leaflet)**
- Phosphatidylethanolamine (POPE): 30%
- Phosphatidylserine (POPS): 15%
- Phosphatidylcholine (POPC): 10%
- Cholesterol (CHL1): 45%

### 1.2 Glycocalyx 구성 요소

RBC glycocalyx는 두께 5-10 nm의 당 사슬 층으로 구성됩니다:

**A. Glycolipids (당지질) - 5-15% of outer leaflet**
- **GM1** (monosialoganglioside): 시알산 1개 함유
- **GM3** (monosialoganglioside): 시알산 1개 함유
- **Globoside (Gb4)**: 중성 당지질

**B. Glycoproteins (당단백질)**
- **Glycophorin A (GPA)**: 
  - RBC 표면의 주요 당단백질 (1백만 개/세포)
  - 16개의 O-glycan 부착
  - 시알산이 풍부 (음전하 제공)
  
- **Band 3**: 
  - 막관통 단백질
  - N-glycan 부착
  
- **CD44**: 
  - Hyaluronan 수용체
  - N-glycan 부착

**C. 시알산 (Sialic Acid)**
- 음전하 (-1) 제공
- Tripod의 L-glucose와 정전기적 반발 유도
- RBC 표면 전하 밀도: -0.02 to -0.05 C/m²

---

## 🛠️ Phase 2: CHARMM-GUI 모델 생성 전략

### Strategy A: Glycolipid-Rich Membrane (추천 - 간단)

**장점**: 
- CHARMM-GUI Membrane Builder에서 직접 생성 가능
- 시알산 함유 당지질로 음전하 표면 구현
- 빠른 시뮬레이션 (glycoprotein보다 작음)

**단점**: 
- Glycoprotein의 입체 장벽 효과 미반영
- 실제 RBC보다 glycocalyx 두께가 얇음

**CHARMM-GUI 설정**:

1. **Membrane Builder** → **Membrane Only System**
2. **Lipid Composition** 설정:
   ```
   Upper Leaflet (Outer):
   - PSM: 25%
   - POPC: 20%
   - GM1: 8%    ← Glycolipid (시알산 함유)
   - GM3: 5%    ← Glycolipid (시알산 함유)
   - CHL1: 42%
   
   Lower Leaflet (Inner):
   - POPE: 30%
   - POPS: 15%
   - POPC: 10%
   - CHL1: 45%
   ```

3. **Box Size**: 
   - X, Y: 120-150 Å (Tripod가 회전할 공간 확보)
   - Z: 충분한 수용액 층 (각 면 30 Å 이상)

4. **Ion Concentration**: 
   - 0.15 M KCl (생리적 조건)

5. **Output**: 
   - OpenMM 또는 GROMACS 선택

---

### Strategy B: Glycoprotein + Glycolipid (현실적 - 복잡)

**장점**: 
- 실제 RBC glycocalyx에 가까움
- 입체 장벽 효과 포함
- 더 두꺼운 glycocalyx 층 (5-10 nm)

**단점**: 
- 복잡한 설정 필요
- 시뮬레이션 시간 증가
- Glycoprotein 구조 준비 필요

**단계별 프로세스**:

#### Step 1: Glycoprotein 준비

```bash
# Glycophorin A 구조 다운로드 (예시)
# PDB에서 glycosylated protein 검색 또는
# CHARMM-GUI Glycan Reader 사용
```

**Option 1**: 기존 glycoprotein PDB 사용
- PDB에서 glycosylated transmembrane protein 검색
- 또는 AlphaFold로 Glycophorin A 구조 생성 후 glycan 추가

**Option 2**: CHARMM-GUI Glycan Reader 사용
1. 단백질 구조 업로드
2. Glycosylation site 지정 (Asn, Ser, Thr)
3. Glycan 구조 선택:
   - **O-glycan**: Core 1 (Galβ1-3GalNAc) + Sialic acid
   - **N-glycan**: Complex type + Sialic acid

#### Step 2: Membrane Builder에 통합

1. **Protein/Membrane System** 선택
2. Glycoprotein 업로드
3. Membrane 조성 설정 (Strategy A와 동일)
4. Glycoprotein 배치:
   - 5-10개의 Glycophorin A 분자 배치
   - 균일하게 분산

---

### Strategy C: Coarse-Grained Model (빠른 스크리닝)

**장점**: 
- 매우 빠른 시뮬레이션 (100x 속도)
- 큰 시스템 가능 (수백 nm)
- Long-timescale dynamics 관찰

**단점**: 
- 원자 수준 상호작용 손실
- Binding energy 정확도 낮음

**사용 도구**: 
- MARTINI force field
- CHARMM-GUI Martini Maker

---

## 🧪 Phase 3: 시뮬레이션 프로토콜

### 3.1 시스템 준비

**A. RBC 모델 (Glycocalyx 있음)**
```python
system_rbc = {
    "membrane": "GM1/GM3-rich (Strategy A)",
    "protein": "GLUT1 (4PYP)",
    "ligand": "Tripod (L-glucose-PEG2)3-TRIS",
    "box_size": [120, 120, 150],  # Å
    "water": "TIP3P",
    "ions": "0.15 M KCl"
}
```

**B. 종양세포 모델 (Glycocalyx 약함)**
```python
system_tumor = {
    "membrane": "Standard POPC/POPE/CHL1 (no glycolipids)",
    "protein": "GLUT1 (4PYP)",
    "ligand": "Tripod",
    "box_size": [120, 120, 150],
    "water": "TIP3P",
    "ions": "0.15 M KCl"
}
```

**C. 대조군 (Naked)**
```python
system_control = {
    "membrane": "POPC only (simple)",
    "protein": "GLUT1 (4PYP)",
    "ligand": "Tripod",
    "box_size": [100, 100, 120],
    "water": "TIP3P",
    "ions": "0.15 M KCl"
}
```

### 3.2 시뮬레이션 단계

#### Stage 1: Equilibration (평형화)
```bash
# 1. Energy Minimization
# 2. NVT equilibration (100 ps, 310 K)
# 3. NPT equilibration (1 ns, 1 bar)
# 4. Production-ready check
```

#### Stage 2: Ligand Approach (접근 시뮬레이션)
```python
# Tripod를 막 위 20-30 Å에 배치
# 자유 확산 시뮬레이션 (50-100 ns)
# 
# 측정 항목:
# - Tripod와 막 표면 간 거리
# - Glycocalyx 침투 여부
# - GLUT1 입구 도달 시간
```

#### Stage 3: Binding Simulation (결합 시뮬레이션)
```python
# Tripod를 GLUT1 입구 근처에 배치
# Production run (100 ns x 3 replicates)
#
# 측정 항목:
# - RMSD (ligand position)
# - Contact analysis (Tripod-GLUT1)
# - Hydrogen bonds
# - MM/PBSA binding energy
```

### 3.3 분석 지표

**A. Glycocalyx 장벽 효과**
```python
metrics = {
    "penetration_depth": "Tripod가 막에 얼마나 가까이 접근했는가?",
    "contact_time": "GLUT1과 접촉한 시간 비율",
    "binding_events": "결합 이벤트 발생 횟수",
    "residence_time": "GLUT1에 머문 평균 시간"
}
```

**예상 결과**:
- **RBC (Glycocalyx 있음)**: 
  - Penetration depth: > 15 Å (막에 도달 못함)
  - Contact time: < 5%
  - Binding events: 0-1
  
- **종양세포 (Glycocalyx 약함)**: 
  - Penetration depth: < 5 Å (막 표면 도달)
  - Contact time: 30-50%
  - Binding events: 3-5

**B. 정전기적 반발**
```python
# Tripod (중성/약양성) vs Glycocalyx (음전하)
# Electrostatic potential map 생성
# PMF (Potential of Mean Force) 계산
```

**C. 입체 장벽**
```python
# Glycan 사슬의 공간 점유율
# Tripod의 회전 반경과 glycocalyx 밀도 비교
```

---

## 📝 Phase 4: Python 스크립트 준비

### 4.1 시스템 설정 스크립트

```python
# scripts/setup_glycocalyx_system.py
"""
CHARMM-GUI 출력 파일을 OpenMM으로 변환하고
Tripod ligand를 추가하는 스크립트
"""

import openmm as mm
import openmm.app as app
from openmm import unit
import mdtraj as md
import numpy as np

def load_charmm_system(psf_file, pdb_file, toppar_dir):
    """CHARMM-GUI 출력 로드"""
    psf = app.CharmmPsfFile(psf_file)
    pdb = app.PDBFile(pdb_file)
    
    # Force field 로드
    params = app.CharmmParameterSet(
        f'{toppar_dir}/par_all36_prot.prm',
        f'{toppar_dir}/par_all36_lipid.prm',
        f'{toppar_dir}/par_all36_carb.prm',
        f'{toppar_dir}/toppar_water_ions.str'
    )
    
    system = psf.createSystem(
        params,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.2*unit.nanometer,
        constraints=app.HBonds
    )
    
    return system, psf, pdb

def add_tripod_ligand(pdb, ligand_pdb, position):
    """Tripod ligand를 시스템에 추가"""
    # Ligand 로드
    ligand = md.load(ligand_pdb)
    
    # 위치 설정 (막 위 20 Å)
    ligand.xyz += position
    
    # 시스템에 병합
    combined = pdb.join(ligand)
    
    return combined

def setup_simulation(system, pdb, temperature=310):
    """시뮬레이션 설정"""
    integrator = mm.LangevinMiddleIntegrator(
        temperature*unit.kelvin,
        1.0/unit.picosecond,
        2.0*unit.femtosecond
    )
    
    simulation = app.Simulation(
        pdb.topology,
        system,
        integrator
    )
    
    simulation.context.setPositions(pdb.positions)
    
    return simulation

if __name__ == "__main__":
    # RBC glycocalyx 시스템 로드
    system, psf, pdb = load_charmm_system(
        'step5_assembly.psf',
        'step5_assembly.pdb',
        'toppar'
    )
    
    # Tripod 추가
    combined_pdb = add_tripod_ligand(
        pdb,
        'tripod_peg2_l_glucose.pdb',
        position=[0, 0, 20]  # 막 위 20 Å
    )
    
    # 시뮬레이션 준비
    simulation = setup_simulation(system, combined_pdb)
    
    print("System ready for equilibration")
```

### 4.2 분석 스크립트

```python
# scripts/analyze_glycocalyx_barrier.py
"""
Glycocalyx 장벽 효과 분석
"""

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

def calculate_penetration_depth(traj, ligand_selection, membrane_selection):
    """
    Tripod가 막에 얼마나 가까이 접근했는지 계산
    """
    ligand = traj.atom_slice(traj.top.select(ligand_selection))
    membrane = traj.atom_slice(traj.top.select(membrane_selection))
    
    # Z 좌표 (막에 수직)
    ligand_z = ligand.xyz[:, :, 2].mean(axis=1)
    membrane_z = membrane.xyz[:, :, 2].mean(axis=1)
    
    # 거리 계산
    distance = ligand_z - membrane_z
    
    return distance

def calculate_contact_time(traj, ligand_selection, protein_selection, cutoff=0.5):
    """
    Tripod가 GLUT1과 접촉한 시간 비율
    """
    contacts = md.compute_contacts(
        traj,
        contacts='all',
        scheme='closest-heavy'
    )
    
    contact_frames = np.sum(contacts[0] < cutoff, axis=1) > 0
    contact_ratio = np.sum(contact_frames) / len(traj)
    
    return contact_ratio

def plot_penetration_profile(distances, title):
    """
    침투 깊이 시계열 플롯
    """
    plt.figure(figsize=(10, 6))
    plt.plot(distances, linewidth=1)
    plt.axhline(y=15, color='r', linestyle='--', label='Glycocalyx barrier')
    plt.axhline(y=5, color='g', linestyle='--', label='Membrane surface')
    plt.xlabel('Frame')
    plt.ylabel('Distance from membrane (Å)')
    plt.title(title)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(f'{title.replace(" ", "_")}.png', dpi=300)
    plt.close()

def compare_systems(traj_rbc, traj_tumor, traj_control):
    """
    세 시스템 비교
    """
    results = {}
    
    for name, traj in [('RBC', traj_rbc), ('Tumor', traj_tumor), ('Control', traj_control)]:
        # 침투 깊이
        distance = calculate_penetration_depth(
            traj,
            'resname LIG',  # Tripod
            'resname POPC or resname GM1 or resname GM3'
        )
        
        # 접촉 시간
        contact = calculate_contact_time(
            traj,
            'resname LIG',
            'protein and resname GLU1'
        )
        
        results[name] = {
            'mean_distance': distance.mean(),
            'min_distance': distance.min(),
            'contact_ratio': contact
        }
        
        # 플롯
        plot_penetration_profile(distance, f'{name} System')
    
    return results

if __name__ == "__main__":
    # 트라젝토리 로드
    traj_rbc = md.load('rbc_production.dcd', top='rbc_system.pdb')
    traj_tumor = md.load('tumor_production.dcd', top='tumor_system.pdb')
    traj_control = md.load('control_production.dcd', top='control_system.pdb')
    
    # 분석
    results = compare_systems(traj_rbc, traj_tumor, traj_control)
    
    # 결과 출력
    print("\n=== Glycocalyx Barrier Effect ===")
    for system, metrics in results.items():
        print(f"\n{system}:")
        print(f"  Mean distance: {metrics['mean_distance']:.2f} Å")
        print(f"  Min distance: {metrics['min_distance']:.2f} Å")
        print(f"  Contact ratio: {metrics['contact_ratio']:.2%}")
```

---

## 🎯 Phase 5: 예상 결과 및 해석

### 5.1 성공 기준

**Hypothesis Validation**:

| 시스템 | 침투 깊이 | 접촉 시간 | 결합 에너지 | 해석 |
|--------|-----------|-----------|-------------|------|
| **RBC (Glycocalyx)** | > 15 Å | < 10% | N/A | ✅ 장벽 효과 확인 |
| **종양세포** | < 5 Å | > 40% | -30 kcal/mol | ✅ 선택적 결합 |
| **대조군** | < 3 Å | > 60% | -35 kcal/mol | ✅ 기본 결합력 |

### 5.2 추가 검증

**A. Glycocalyx 밀도 변화**
- GM1/GM3 비율을 5%, 10%, 15%로 변화
- 장벽 효과와 당지질 밀도의 상관관계 확인

**B. Tripod 변형 테스트**
- PEG 길이 변화 (PEG2 vs PEG4 vs PEG6)
- 긴 PEG가 glycocalyx를 더 잘 침투하는지 확인

**C. 경쟁 실험**
- D-glucose 50개 추가
- Tripod vs D-glucose 경쟁 시뮬레이션

---

## 📚 참고 자료

### CHARMM-GUI Tutorials
- Membrane Builder: http://www.charmm-gui.org/?doc=tutorial&project=membrane
- Glycan Reader: http://www.charmm-gui.org/?doc=tutorial&project=glycan

### 논문
1. **RBC Glycocalyx**:
   - Bäumler et al. (2001) "Electrophoresis of human red blood cells"
   - Pries et al. (2000) "The endothelial surface layer"

2. **Ganglioside Structure**:
   - Schnaar et al. (2014) "Glycosphingolipids in neural function"
   - Yu et al. (2012) "Gangliosides in cancer"

3. **Tumor Glycocalyx**:
   - Paszek et al. (2014) "The cancer glycocalyx"
   - Woods et al. (2013) "Tumor glycocalyx and metastasis"

---

## ✅ Next Steps

1. **CHARMM-GUI에서 RBC 모델 생성**
   - Strategy A (Glycolipid-rich) 사용
   - GM1 8%, GM3 5% 포함
   
2. **종양세포 모델 생성**
   - Standard lipid composition
   - No glycolipids
   
3. **시뮬레이션 실행**
   - 각 시스템 100 ns x 3 replicates
   
4. **결과 분석**
   - Penetration depth 비교
   - Contact time 비교
   - 통계적 유의성 검증

---

**질문이나 도움이 필요하면 언제든 알려주세요!**
