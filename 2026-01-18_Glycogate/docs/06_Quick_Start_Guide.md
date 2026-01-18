# OpenMM MD 시뮬레이션 빠른 시작 가이드

**작성일:** 2026-01-18  
**목적:** 1-Arm PEG24-Glc OpenMM MD 빠른 실행

---

## ✅ 준비 완료된 것

1. ✅ 1-arm PEG24-Glc 구조 생성
2. ✅ CHARMM-GUI Ligand Reader (topology)
3. ✅ CHARMM-GUI Solution Builder (수용액 시스템)
4. ✅ OpenMM 파일 생성

**시스템 정보:**
- 리간드: 1-arm TRIS-PEG24-Glc (309 원자)
- 물: TIP3P (~40,000개)
- 이온: 150 mM NaCl
- Box: 12 × 12 × 12 nm³ (120 Å)
- 총 원자 수: ~120,000개

---

## 🚀 MD 시뮬레이션 실행

### 방법 1: 간단 스크립트 (권장)

```bash
cd /home/pjho3/projects/Drug/2026-01-18_Glycogate

# Drug-MD 환경 활성화
conda activate Drug-MD

# Replica 1 실행 (테스트: 10 ns)
bash scripts/09_run_openmm_simple.sh 1

# Replica 2, 3 (병렬 실행 가능)
bash scripts/09_run_openmm_simple.sh 2 &
bash scripts/09_run_openmm_simple.sh 3 &
```

**실행 내용:**
1. 에너지 최소화 (5000 steps)
2. NVT 평형화 (125 ps, 1 fs timestep)
3. Production MD (10 ns, 4 fs timestep) - 테스트용

**출력 파일:**
```
results/md_1arm_openmm/
├── step4_equilibration_rep1.rst    # 평형화 restart
├── step4_equilibration_rep1.dcd    # 평형화 trajectory
├── step4_equilibration_rep1.out    # 평형화 로그
├── step5_1_rep1.rst                # Production step 1
├── step5_1_rep1.dcd                # Production trajectory 1
├── step5_1_rep1.out                # Production 로그 1
├── step5_2_rep1.rst                # Production step 2
└── ...
```

---

### 방법 2: 수동 실행 (고급)

```bash
cd /home/pjho3/projects/Drug/2026-01-18_Glycogate/data/solution\ builder/openmm

# 1. 평형화
python -u openmm_run.py \
    -i step4_equilibration.inp \
    -t toppar.str \
    -p step3_input.psf \
    -c step3_input.crd \
    -b sysinfo.dat \
    -orst equilibration.rst \
    -odcd equilibration.dcd \
    > equilibration.out

# 2. Production MD (1 ns)
python -u openmm_run.py \
    -i step5_production.inp \
    -t toppar.str \
    -p step3_input.psf \
    -c step3_input.crd \
    -irst equilibration.rst \
    -orst production_1.rst \
    -odcd production_1.dcd \
    > production_1.out
```

---

## 📊 시뮬레이션 파라미터

### 평형화 (step4_equilibration.inp)

```
에너지 최소화: 5000 steps
NVT 평형화: 125 ps (125,000 steps × 1 fs)
온도: 303.15 K (30°C)
압력: 없음 (NVT)
Restraints: Yes (backbone 400, sidechain 40 kJ/mol/nm²)
```

### Production (step5_production.inp)

```
시간: 1 ns per step (250,000 steps × 4 fs)
온도: 303.15 K
압력: 1 bar (NPT, isotropic)
Restraints: No
DCD 저장: 100 ps 간격
```

**200 ns 시뮬레이션:**
- 200번 반복 필요
- 스크립트에서 `CNTMAX=200`으로 변경

---

## 🔧 200 ns 실행으로 변경

`scripts/09_run_openmm_simple.sh` 파일 수정:

```bash
# 라인 73-74 수정
# CNTMAX=10  # 10 ns (테스트용)
CNTMAX=200  # 200 ns (실제 실행용)
```

그 다음 실행:
```bash
bash scripts/09_run_openmm_simple.sh 1
```

**예상 시간 (GPU 기준):**
- NVIDIA RTX 3090: ~1.5-2일
- NVIDIA A100: ~1일

---

## 📈 실시간 모니터링

### 로그 확인

```bash
# 평형화 로그
tail -f results/md_1arm_openmm/step4_equilibration_rep1.out

# Production 로그
tail -f results/md_1arm_openmm/step5_1_rep1.out
```

### 진행률 확인

```bash
# 완료된 스텝 확인
ls results/md_1arm_openmm/step5_*_rep1.rst | wc -l

# 총 200 스텝 중 몇 개 완료되었는지 확인
```

### 에너지 플롯 (Python)

```python
import re
import matplotlib.pyplot as plt

# 로그 파일 읽기
with open('results/md_1arm_openmm/step5_1_rep1.out') as f:
    lines = f.readlines()

# 에너지 추출
times = []
energies = []
for line in lines:
    if 'Progress' in line:
        # 예: Progress: 10.0%, Time: 0.1 ps, Speed: 150 ns/day
        match = re.search(r'Time: ([\d.]+)', line)
        if match:
            times.append(float(match.group(1)))

# 플롯
plt.plot(times)
plt.xlabel('Step')
plt.ylabel('Time (ps)')
plt.savefig('progress.png')
```

---

## 🔍 문제 해결

### 문제 1: OpenMM 없음

```bash
conda activate Drug-MD
conda install -c conda-forge openmm
```

### 문제 2: CUDA 없음

**증상:**
```
Platform 'CUDA' not found
```

**해결:**
- CPU로 실행됨 (매우 느림)
- GPU 서버 사용 권장

### 문제 3: 메모리 부족

**증상:**
```
CUDA out of memory
```

**해결:**
- `openmm_run.py`에서 precision 변경
- 또는 더 큰 GPU 사용

### 문제 4: 시뮬레이션 불안정

**증상:**
- NaN 에너지
- 구조 폭발

**해결:**
1. 평형화 더 길게 (step4_equilibration.inp의 nstep 증가)
2. Timestep 줄이기 (dt = 0.002)
3. 초기 구조 확인

---

## 📊 다음 단계: 분석

시뮬레이션 완료 후:

1. **Trajectory 병합**
   ```bash
   # MDTraj 사용
   python scripts/10_merge_trajectories.py
   ```

2. **End-to-end 거리 분석**
   ```python
   import mdtraj as md
   traj = md.load('merged.dcd', top='step3_input.pdb')
   # 분석 코드...
   ```

3. **Radius of gyration**
4. **Tail 확률 계산**
5. **3 replica 통합 분석**

---

## 📝 체크리스트

실행 전:
- [ ] Drug-MD 환경 활성화
- [ ] OpenMM 설치 확인
- [ ] GPU 사용 가능 확인
- [ ] 디스크 공간 확인 (~50 GB per replica)

실행 중:
- [ ] 로그 파일 모니터링
- [ ] 에너지 안정성 확인
- [ ] 진행률 추적

실행 후:
- [ ] Trajectory 파일 확인
- [ ] 분석 스크립트 실행
- [ ] 결과 판정

---

**작성자:** Cascade AI  
**최종 수정:** 2026-01-18
