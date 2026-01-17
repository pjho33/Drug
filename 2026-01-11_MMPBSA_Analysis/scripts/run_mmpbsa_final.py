#!/usr/bin/env python3
"""
MMPBSA.py 최종 실행
"""

import subprocess
import os
from pathlib import Path

# 출력 디렉토리 읽기
with open('/tmp/mmpbsa_outdir.txt', 'r') as f:
    OUTDIR = Path(f.read().strip())

AMBER_DIR = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/GLUT1SDGComplex260110/amber")

print("=" * 80)
print("MMPBSA.py 실행")
print("=" * 80)
print()
print(f"작업 디렉토리: {OUTDIR}")
print(f"좌표 파일: {AMBER_DIR / 'step5_input.rst7'}")
print()

# MMPBSA 입력 파일 생성
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

input_file = OUTDIR / "mmpbsa.in"
with open(input_file, 'w') as f:
    f.write(mmpbsa_input)

print("✅ MMPBSA 입력 파일 생성")
print()

# MMPBSA.py 실행
print("=" * 80)
print("MMPBSA.py 실행 중...")
print("=" * 80)
print()

RADII_PARM = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/mmpbsa_amber/step5_input_radii.parm7")

cmd = [
    "MMPBSA.py",
    "-O",
    "-i", str(input_file),
    "-o", str(OUTDIR / "FINAL_RESULTS_MMPBSA.dat"),
    "-sp", str(RADII_PARM),  # 원본 radii 포함 topology (140,250 atoms)
    "-cp", str(OUTDIR / "complex.prmtop"),
    "-rp", str(OUTDIR / "receptor.prmtop"),
    "-lp", str(OUTDIR / "ligand.prmtop"),
    "-y", str(AMBER_DIR / "step5_input.rst7")
]

print("명령:")
print(" ".join(cmd))
print()

result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(OUTDIR))

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
    
    results_file = OUTDIR / "FINAL_RESULTS_MMPBSA.dat"
    if results_file.exists():
        print("결과:")
        print("-" * 80)
        with open(results_file, 'r') as f:
            print(f.read())
    else:
        print("⚠️  결과 파일이 생성되지 않았습니다.")
else:
    print(f"❌ MMPBSA.py 실패 (exit code: {result.returncode})")
    print("=" * 80)
    exit(1)
