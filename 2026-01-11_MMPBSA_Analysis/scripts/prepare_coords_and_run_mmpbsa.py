#!/usr/bin/env python3
"""
좌표 파일 준비 및 MMPBSA.py 실행
"""

import subprocess
import os
from pathlib import Path

# 출력 디렉토리 (parmed_out)
OUTDIR = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/mmpbsa_amber/parmed_out")
AMBER_DIR = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/GLUT1SDGComplex260110/amber")
RADII_PARM = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/mmpbsa_amber/step5_input_radii.parm7")

os.chdir(OUTDIR)

print("=" * 80)
print("좌표 파일 준비 및 MMPBSA.py 실행")
print("=" * 80)
print()
print(f"작업 디렉토리: {OUTDIR}")
print()

# Step 1: cpptraj로 물/이온 제거한 좌표 생성
print("Step 1: 물/이온 제거한 좌표 생성")
print("-" * 80)

cpptraj_input = f"""parm {RADII_PARM}
trajin {AMBER_DIR / 'step5_input.rst7'}
strip :POT,CLA,WAT
trajout complex.rst7 restart
"""

with open("strip_solvent.in", "w") as f:
    f.write(cpptraj_input)

result = subprocess.run(
    ["cpptraj", "-i", "strip_solvent.in"],
    capture_output=True,
    text=True
)

print(result.stdout)
if result.returncode != 0:
    print("STDERR:", result.stderr)
    print("❌ cpptraj 실패")
    exit(1)

print("✅ 좌표 파일 생성 완료: complex.rst7")
print()

# Step 2: MMPBSA 입력 파일 생성
print("Step 2: MMPBSA 입력 파일 생성")
print("-" * 80)

mmpbsa_input = """MMPBSA calculation with PB
&general
startframe=1
endframe=1
interval=1
verbose=2
/

&pb
istrng=0.150
fillratio=4.0
/
"""

with open("mmpbsa.in", "w") as f:
    f.write(mmpbsa_input)

print("✅ mmpbsa.in 생성 완료")
print()

# Step 3: MMPBSA.py 실행
print("Step 3: MMPBSA.py 실행")
print("=" * 80)

cmd = [
    "MMPBSA.py",
    "-O",
    "-i", "mmpbsa.in",
    "-o", "FINAL_RESULTS_MMPBSA.dat",
    "-sp", "complex.prmtop",  # 물/이온 제거된 complex
    "-cp", "complex.prmtop",
    "-rp", "receptor.prmtop",
    "-lp", "ligand.prmtop",
    "-y", "complex.rst7"
]

print("명령:")
print(" ".join(cmd))
print()

result = subprocess.run(cmd, capture_output=True, text=True)

print(result.stdout)
if result.stderr:
    print("STDERR:")
    print(result.stderr)

print()
print("=" * 80)
if result.returncode == 0:
    print("✅ MMPBSA.py 성공!")
    print("=" * 80)
    print()
    
    if Path("FINAL_RESULTS_MMPBSA.dat").exists():
        print("결과:")
        print("-" * 80)
        with open("FINAL_RESULTS_MMPBSA.dat", "r") as f:
            print(f.read())
    else:
        print("⚠️  결과 파일이 생성되지 않았습니다.")
else:
    print(f"❌ MMPBSA.py 실패 (exit code: {result.returncode})")
    print("=" * 80)
    exit(1)
