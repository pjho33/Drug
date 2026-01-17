# CHARMM-GUI Membrane Builder 제출 가이드

## 📁 준비된 파일

### 실험군 (Glycosylated GLUT1 + Tripod)
```
/home/pjho3tr/projects/Drug/phase3_charmm_gui_submission/experimental/
├── glut1_tripod_complex.pdb       (7,702 atoms)
├── trp.rtf                         (Tripod topology)
└── trp.prm                         (Tripod parameters)
```

### 대조군 (Non-glycosylated GLUT1 + Tripod)
```
/home/pjho3tr/projects/Drug/phase3_charmm_gui_submission/control/
├── glut1_tripod_complex_control.pdb  (7,114 atoms)
├── trp.rtf                            (Tripod topology)
└── trp.prm                            (Tripod parameters)
```

---

## 🌐 제출 절차

### 1️⃣ CHARMM-GUI 접속
**URL**: https://www.charmm-gui.org/?doc=input/membrane

또는:
1. https://www.charmm-gui.org 접속
2. **Input Generator** 클릭
3. **Membrane Builder** 선택

---

### 2️⃣ 실험군 제출 (먼저)

#### Step 1: PDB Upload
- **Upload PDB**: `experimental/glut1_tripod_complex.pdb` 선택
- **Next Step** 클릭

#### Step 2: Manipulate PDB
- 특별한 수정 없이 **Next Step**

#### Step 3: Identify Components
- **Protein**: 자동 인식됨
- **Ligand/Cofactor**: TRP (Tripod) 확인
  - **Topology file**: `experimental/trp.rtf` 업로드
  - **Parameter file**: `experimental/trp.prm` 업로드
- **Glycans**: 자동 인식 확인 (Asn45, Asn247, Asn255)
- **Next Step**

#### Step 4: Orient Molecule
- **Orientation**: Membrane Builder가 자동 계산
- 또는 기존 설정 참고:
  - Z-axis: Membrane normal
- **Next Step**

#### Step 5: Determine System Size
- **System Size**:
  ```
  X: 100 Å (또는 자동 권장값)
  Y: 100 Å
  Z: 100 Å
  ```
- **Membrane type**: POPC (또는 POPC:POPE 7:3)
- **Next Step**

#### Step 6: Build Components
- **Solvate**: Yes
- **Water model**: TIP3P
- **Ion concentration**: 0.15 M KCl
- **Neutralize system**: Yes
- **Next Step**

#### Step 7: Equilibration
- **Skip equilibration** (우리가 직접 할 것)
- **Next Step**

#### Step 8: Output
- **Force field**: CHARMM36m
- **Output format**: **OpenMM** ✅
- **Job name**: `glut1_tripod_glycosylated`
- **Submit**

---

### 3️⃣ 대조군 제출 (동일 과정)

위 과정을 `control/glut1_tripod_complex_control.pdb`로 반복

**차이점**:
- Glycans 없음 (자동 인식됨)
- Job name: `glut1_tripod_control`

---

## ⏱️ 대기 시간

- **예상**: 2-3시간
- **이메일 알림**: 완료 시 자동 발송
- **다운로드**: 링크 클릭하여 ZIP 파일 다운로드

---

## 📥 다운로드 후 실행

### 압축 해제
```bash
cd ~/Downloads
unzip charmm-gui-*glycosylated*.zip -d ~/Downloads/experimental_tripod
unzip charmm-gui-*control*.zip -d ~/Downloads/control_tripod
```

### 즉시 실행 (GPU 0, 1 동시)
```bash
cd ~/projects/Drug/phase3_with_tripod

# 실험군 스크립트 생성
cat > run_exp_final.py << 'PYEOF'
import sys, os
sys.path.insert(0, os.path.expanduser('~/Downloads/experimental_tripod/openmm'))
os.chdir(os.path.expanduser('~/Downloads/experimental_tripod/openmm'))

from omm_run import main
sys.argv = ['openmm_run.py', '--platform', 'CUDA', '--device', '0']
main()
PYEOF

# 대조군 스크립트 생성
cat > run_ctrl_final.py << 'PYEOF'
import sys, os
sys.path.insert(0, os.path.expanduser('~/Downloads/control_tripod/openmm'))
os.chdir(os.path.expanduser('~/Downloads/control_tripod/openmm'))

from omm_run import main
sys.argv = ['openmm_run.py', '--platform', 'CUDA', '--device', '1']
main()
PYEOF

# 동시 실행
python run_exp_final.py > exp_final.log 2>&1 &
python run_ctrl_final.py > ctrl_final.log 2>&1 &

echo "Simulations started!"
```

---

## ✅ 체크리스트

제출 전 확인:
- [ ] 실험군 PDB 파일 준비됨
- [ ] 대조군 PDB 파일 준비됨
- [ ] trp.rtf, trp.prm 파일 준비됨
- [ ] CHARMM-GUI 계정 로그인됨

제출 시:
- [ ] Ligand topology/parameter 업로드 확인
- [ ] Glycans 인식 확인 (실험군만)
- [ ] Output format: OpenMM 선택
- [ ] Job name 구분 가능하게 설정

---

## 🎯 예상 결과

다운로드 파일 구조:
```
charmm-gui-*/
├── openmm/
│   ├── step5_input.psf      ← Tripod 포함 topology
│   ├── step5_input.pdb      ← Tripod 포함 좌표
│   ├── step5_input.crd
│   ├── toppar.str           ← 모든 parameters 포함
│   ├── openmm_run.py        ← 실행 스크립트
│   └── ...
└── toppar/                  ← 모든 force field 파일
```

**즉시 실행 가능!** 추가 수정 불필요!

---

## 📞 문제 발생 시

1. **Ligand 인식 안 됨**: trp.rtf에서 residue name이 TRP인지 확인
2. **Parameter 오류**: trp.prm 형식 확인
3. **제출 실패**: 파일 크기 확인 (10MB 이하)

---

**지금 바로 제출하세요!**
https://www.charmm-gui.org/?doc=input/membrane
