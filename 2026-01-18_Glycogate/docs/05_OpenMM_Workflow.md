# OpenMM 워크플로우

**작성일:** 2026-01-18  
**목적:** 1-Arm PEG24-Glc OpenMM MD 시뮬레이션

---

## 🎯 OpenMM vs GROMACS

### OpenMM의 장점

| 항목 | OpenMM | GROMACS |
|------|--------|---------|
| **언어** | Python | C++ (명령행) |
| **스크립트** | 쉬움 | 복잡 (bash) |
| **GPU 성능** | 우수 | 우수 |
| **유연성** | 매우 높음 | 제한적 |
| **분석 통합** | 쉬움 (MDTraj) | 별도 도구 |
| **커스터마이징** | Python 코드 | 소스 수정 |
| **CHARMM-GUI** | 지원 | 지원 |

**결론:** Python 기반 워크플로우에 OpenMM이 더 적합

---

## 📋 전체 워크플로우

### 1. CHARMM-GUI Solution Builder (OpenMM 출력)

**URL:** http://www.charmm-gui.org/?doc=input/solution

**절차:**
1. **Upload PDB**
   - 파일: `data/charmm-gui-6871763698/ligandrm.pdb`

2. **System Setup**
   - Water model: TIP3P
   - Ion concentration: 150 mM NaCl
   - Box type: Cubic
   - Box size: 12.0 nm
   - Neutralize: Yes

3. **Output Selection**
   - **OpenMM** 선택 ⭐
   - Download

4. **파일 구성**
   ```
   charmm-gui-XXXXXX/openmm/
   ├── step3_input.pdb       # 초기 구조 (물+이온 포함)
   ├── step3_input.psf       # PSF topology
   ├── toppar/               # CHARMM36 force field
   │   ├── par_all36m_prot.prm
   │   ├── par_all36_cgenff.prm
   │   └── ...
   ├── openmm_run.py         # CHARMM-GUI 제공 스크립트
   └── step5_*.inp           # 설정 파일
   ```

---

### 2. OpenMM 환경 설정

**Drug-MD 환경 확인:**
```bash
conda activate Drug-MD
python -c "import openmm; print(openmm.version.version)"
```

**설치 (필요시):**
```bash
conda install -c conda-forge openmm
conda install -c conda-forge mdtraj
conda install -c conda-forge mdanalysis
```

**GPU 확인:**
```python
import openmm
print(openmm.Platform.getPlatformByName('CUDA'))
```

---

### 3. MD 시뮬레이션 실행

#### 방법 1: 제공된 스크립트 사용 (권장)

```bash
cd /home/pjho3/projects/Drug/2026-01-18_Glycogate
conda activate Drug-MD

# 스크립트 수정 (CHARMM-GUI ID 업데이트)
# scripts/08_run_openmm_md.py의 charmm_gui_dir 경로 수정

# Replica 1 실행
python scripts/08_run_openmm_md.py 1

# Replica 2, 3 (병렬 실행 가능)
python scripts/08_run_openmm_md.py 2 &
python scripts/08_run_openmm_md.py 3 &
```

**출력 파일:**
```
results/md_1arm_openmm/
├── md_rep1.dcd              # Trajectory (10 ps 간격)
├── md_rep1.log              # 에너지, 온도 등
├── md_rep1.chk              # Checkpoint (재시작용)
├── md_rep1_equilibrated.pdb # 평형화 후 구조
├── md_rep2.*
└── md_rep3.*
```

---

#### 방법 2: CHARMM-GUI 제공 스크립트 수정

CHARMM-GUI의 `openmm_run.py`를 수정하여 사용:

```python
# openmm_run.py 수정 예시

# 시뮬레이션 길이 변경
nsteps = 100000000  # 200 ns (2 fs timestep)

# Reporter 간격 조정
dcd_reporter = DCDReporter('output.dcd', 5000)  # 10 ps
state_reporter = StateDataReporter('output.log', 5000)

# Replica 별 random seed
integrator.setRandomNumberSeed(replica_number)
```

---

### 4. 시뮬레이션 모니터링

**실시간 모니터링:**
```bash
# 로그 파일 확인
tail -f results/md_1arm_openmm/md_rep1.log

# 진행률 확인
grep "Progress" results/md_1arm_openmm/md_rep1.log | tail -1
```

**에너지 플롯 (Python):**
```python
import pandas as pd
import matplotlib.pyplot as plt

# 로그 파일 읽기
data = pd.read_csv('md_rep1.log', sep=',', skiprows=1)

# 에너지 플롯
plt.figure(figsize=(10, 6))
plt.plot(data['Time (ps)'], data['Potential Energy (kJ/mole)'])
plt.xlabel('Time (ps)')
plt.ylabel('Potential Energy (kJ/mol)')
plt.savefig('energy.png')
```

---

### 5. 궤적 분석

#### MDTraj 사용

```python
import mdtraj as md
import numpy as np

# Trajectory 로드
traj = md.load('md_rep1.dcd', top='md_rep1_equilibrated.pdb')

# End-to-end 거리 (TRIS 중심 - Glucose 말단)
# 원자 인덱스는 구조에 따라 조정 필요
tris_center = traj.topology.select('resname LIG and name C1')  # 예시
glc_end = traj.topology.select('resname LIG and name O6')      # 예시

distances = md.compute_distances(traj, [[tris_center[0], glc_end[0]]])

# Radius of gyration
rg = md.compute_rg(traj)

# 결과 저장
np.save('end_to_end.npy', distances)
np.save('rg.npy', rg)

# 분포 플롯
import matplotlib.pyplot as plt

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.hist(distances * 10, bins=50, density=True)  # nm 단위
plt.xlabel('End-to-end distance (nm)')
plt.ylabel('Probability density')

plt.subplot(1, 2, 2)
plt.hist(rg * 10, bins=50, density=True)
plt.xlabel('Radius of gyration (nm)')
plt.ylabel('Probability density')

plt.tight_layout()
plt.savefig('analysis.png')
```

---

#### MDAnalysis 사용

```python
import MDAnalysis as mda
from MDAnalysis.analysis import rms, distances

# Universe 생성
u = mda.Universe('md_rep1_equilibrated.pdb', 'md_rep1.dcd')

# 리간드 선택
lig = u.select_atoms('resname LIG')

# End-to-end 거리
tris = u.select_atoms('resname LIG and name C1')
glc = u.select_atoms('resname LIG and name O6')

end_to_end = []
for ts in u.trajectory:
    d = distances.dist(tris, glc)[2][0]
    end_to_end.append(d)

# Radius of gyration
rg = []
for ts in u.trajectory:
    rg.append(lig.radius_of_gyration())
```

---

### 6. 분석 스크립트 (자동화)

분석 스크립트 작성 예정:
- `scripts/09_analyze_trajectory.py`
- End-to-end 거리 분포
- Rg 분포
- Tail 확률 계산
- Autocorrelation time
- 3 replica 통합 분석

---

## 🔧 시뮬레이션 파라미터

### 기본 설정

```python
# Integrator
LangevinMiddleIntegrator(
    temperature=300*kelvin,
    frictionCoeff=1.0/picosecond,
    stepSize=0.002*picoseconds  # 2 fs
)

# Barostat (NPT)
MonteCarloBarostat(
    pressure=1.0*bar,
    temperature=300*kelvin,
    frequency=25
)

# Nonbonded
nonbondedMethod=PME
nonbondedCutoff=1.2*nanometer
ewaldErrorTolerance=0.0005

# Constraints
constraints=HBonds
rigidWater=True
```

### 시뮬레이션 길이

| 단계 | 시간 | 스텝 수 |
|------|------|---------|
| 에너지 최소화 | - | 1000 |
| NVT 평형화 | 100 ps | 50,000 |
| Production | 200 ns | 100,000,000 |

**총 3 replica × 200 ns = 600 ns**

---

## 📊 예상 성능

### GPU 성능 (NVIDIA RTX 3090 기준)

- **속도:** ~100-150 ns/day
- **200 ns 시뮬레이션:** 1.5-2일
- **3 replica (병렬):** 1.5-2일

### CPU 성능 (80 코어 기준)

- **속도:** ~5-10 ns/day
- **200 ns 시뮬레이션:** 20-40일
- **권장하지 않음**

---

## 🚨 문제 해결

### 문제 1: CUDA 없음

**증상:**
```
Platform 'CUDA' not found
```

**해결:**
```bash
# CUDA 설치 확인
nvidia-smi

# OpenMM CUDA 재설치
conda install -c conda-forge openmm cudatoolkit=11.8
```

---

### 문제 2: 메모리 부족

**증상:**
```
CUDA out of memory
```

**해결:**
- Precision을 'single'로 변경
- 또는 'mixed' 사용 (기본)
- Reporter 간격 늘리기

---

### 문제 3: 시뮬레이션 불안정

**증상:**
- NaN 에너지
- 폭발하는 구조

**해결:**
1. 에너지 최소화 더 길게
2. NVT 평형화 더 길게 (500 ps)
3. Timestep 줄이기 (1 fs)
4. 초기 구조 확인

---

## 📚 참고 자료

### OpenMM 문서
- http://docs.openmm.org/
- http://docs.openmm.org/latest/userguide/

### CHARMM-GUI
- http://www.charmm-gui.org/
- Solution Builder 튜토리얼

### 분석 도구
- MDTraj: http://mdtraj.org/
- MDAnalysis: https://www.mdanalysis.org/

---

**작성자:** Cascade AI  
**최종 수정:** 2026-01-18
