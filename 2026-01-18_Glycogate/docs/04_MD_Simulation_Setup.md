# 1-Arm MD 시뮬레이션 설정

**작성일:** 2026-01-18  
**목적:** TRIS-PEG24-Glc 1-arm 시스템 MD 시뮬레이션 준비

---

## 🎯 시스템 구성

### 현재 상태
- ✅ CHARMM-GUI로 리간드 topology 생성 완료
- ✅ GROMACS 파일 생성 (`LIG.itp`, `topol.top`)
- ⏳ 수용액 시스템 구축 필요

### 목표 시스템
```
1-arm 리간드 (1개)
    +
TIP3P 물 (~40,000개)
    +
150 mM NaCl
```

**Box size:** 12 × 12 × 12 nm³ (또는 적절한 크기)

---

## 📋 작업 순서

### 1. CHARMM-GUI Solution Builder 사용 (권장)

**URL:** http://www.charmm-gui.org/?doc=input/solution

**절차:**
1. **Upload Structure**
   - PDB 파일 업로드: `ligandrm.pdb`
   - 또는 MOL2: `1arm_peg24_glc.mol2`

2. **System Setup**
   - Water model: TIP3P
   - Ion concentration: 150 mM NaCl
   - Box type: Cubic
   - Box size: 12.0 nm (또는 auto)
   - Neutralize: Yes

3. **Output**
   - GROMACS 선택
   - Download

4. **파일 구성**
   ```
   step3_input.gro    # 초기 구조
   step3_input.pdb    # 초기 구조 (PDB)
   topol.top          # Topology
   *.itp              # Force field
   *.mdp              # MD 파라미터
   ```

---

### 2. 수동 시스템 구축 (GROMACS)

현재 CHARMM-GUI 파일이 리간드만 있으므로 수동으로 물과 이온 추가:

```bash
# 1. Box 생성
gmx editconf -f ligandrm.pdb -o box.gro -c -d 6.0 -bt cubic

# 2. 물 추가
gmx solvate -cp box.gro -cs spc216.gro -o solvated.gro -p topol.top

# 3. 이온 추가 (150 mM NaCl)
gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o system.gro -p topol.top -pname NA -nname CL -conc 0.15 -neutral

# 4. 에너지 최소화
gmx grompp -f em.mdp -c system.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

# 5. NVT 평형화
gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt

# 6. NPT 평형화
gmx grompp -f npt.mdp -c nvt.gro -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt

# 7. Production MD
gmx grompp -f md.mdp -c npt.gro -p topol.top -o md.tpr
gmx mdrun -v -deffnm md
```

---

## 📝 MDP 파일 설정

### Energy Minimization (em.mdp)

```
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000

nstlist     = 10
cutoff-scheme = Verlet
ns_type     = grid
rlist       = 1.2
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
```

### NVT Equilibration (nvt.mdp)

```
integrator  = md
dt          = 0.002
nsteps      = 50000  ; 100 ps

nstxout     = 5000
nstvout     = 5000
nstenergy   = 500
nstlog      = 500

cutoff-scheme = Verlet
nstlist     = 10
rlist       = 1.2
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2

tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.1
ref_t       = 300

pbc         = xyz
gen_vel     = yes
gen_temp    = 300
gen_seed    = -1
```

### NPT Equilibration (npt.mdp)

```
integrator  = md
dt          = 0.002
nsteps      = 50000  ; 100 ps

nstxout     = 5000
nstvout     = 5000
nstenergy   = 500
nstlog      = 500

cutoff-scheme = Verlet
nstlist     = 10
rlist       = 1.2
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2

tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.1
ref_t       = 300

pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5

pbc         = xyz
gen_vel     = no
```

### Production MD (md.mdp)

```
integrator  = md
dt          = 0.002
nsteps      = 100000000  ; 200 ns

nstxout-compressed = 5000  ; 10 ps
compressed-x-grps  = System

nstenergy   = 5000
nstlog      = 5000

cutoff-scheme = Verlet
nstlist     = 10
rlist       = 1.2
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2

tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.1
ref_t       = 300

pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5

pbc         = xyz
gen_vel     = no

; Constraints
constraints     = h-bonds
constraint_algorithm = LINCS
```
1
---

## 🔧 다음 단계

1. **CHARMM-GUI Solution Builder로 시스템 구축** (권장)
   - 또는 수동으로 물/이온 추가

2. **3 Replica 준비**
   - 다른 random seed 사용
   - `gen_seed = -1` (자동), 또는 `42`, `123`, `456`

3. **MD 실행**
   - 200-500 ns × 3 replica
   - GPU 사용 권장

4. **분석**
   - End-to-end 거리
   - Radius of gyration
   - Tail 확률

---

**작성자:** Cascade AI  
**최종 수정:** 2026-01-18
