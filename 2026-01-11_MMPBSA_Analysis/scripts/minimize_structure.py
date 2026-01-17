#!/usr/bin/env python3
"""
구조 Minimization (클래시 제거)
"""

import subprocess
from pathlib import Path
import os

AMBER_DIR = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/GLUT1SDGComplex260110/amber")
RADII_PARM = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/mmpbsa_amber/step5_input_radii.parm7")
RST7 = AMBER_DIR / "step5_input.rst7"
OUTDIR = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/mmpbsa_amber/minimized")

OUTDIR.mkdir(exist_ok=True)
os.chdir(OUTDIR)

print("=" * 80)
print("구조 Minimization (클래시 제거)")
print("=" * 80)
print()

# Minimization 입력 파일 생성
min_input = """Minimization to remove clashes
 &cntrl
  imin=1,           ! Minimization
  maxcyc=5000,      ! Maximum cycles
  ncyc=2500,        ! Steepest descent steps
  ntb=0,            ! No periodic boundary
  igb=8,            ! GB model (implicit solvent)
  cut=999.0,        ! Non-bonded cutoff
  ntpr=100,         ! Print frequency
  ntr=0,            ! No restraints
 /
"""

with open("minimize.in", "w") as f:
    f.write(min_input)

print("Minimization 입력 파일 생성 완료")
print()

# Sander로 minimization 실행
print("=" * 80)
print("Sander Minimization 실행 중...")
print("=" * 80)
print()

cmd = [
    "sander",
    "-O",
    "-i", "minimize.in",
    "-o", "minimize.out",
    "-p", str(RADII_PARM),
    "-c", str(RST7),
    "-r", "minimized.rst7",
    "-ref", str(RST7)
]

print("명령:")
print(" ".join(cmd))
print()

result = subprocess.run(cmd, capture_output=True, text=True)

if result.returncode == 0:
    print("✅ Minimization 성공!")
    print()
    
    # 결과 확인
    if Path("minimize.out").exists():
        with open("minimize.out", "r") as f:
            lines = f.readlines()
            
        # 에너지 변화 확인
        print("에너지 변화:")
        print("-" * 80)
        for line in lines:
            if "NSTEP" in line or "ENERGY" in line or "BOND" in line:
                print(line.strip())
        
        # 최종 에너지 찾기
        for i, line in enumerate(lines):
            if "FINAL RESULTS" in line:
                for j in range(i, min(i+20, len(lines))):
                    print(lines[j].strip())
                break
    
    print()
    print("=" * 80)
    print(f"✅ Minimized 구조: {OUTDIR}/minimized.rst7")
    print("=" * 80)
    
else:
    print("❌ Minimization 실패")
    print()
    print("STDOUT:")
    print(result.stdout)
    print()
    print("STDERR:")
    print(result.stderr)
    exit(1)

# 저장 경로 기록
with open('/tmp/minimized_rst7.txt', 'w') as f:
    f.write(str(OUTDIR / "minimized.rst7"))

print()
print("다음 단계: minimized 구조로 MMPBSA 재계산")
