#!/usr/bin/env python3
"""
GB 방법으로 MMPBSA 재계산 (현재 구조 그대로)
"""

import subprocess
import os
from pathlib import Path

# ParmEd로 생성한 topology 사용
OUTDIR = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/mmpbsa_amber/parmed_out")
AMBER_DIR = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/GLUT1SDGComplex260110/amber")
RADII_PARM = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/mmpbsa_amber/step5_input_radii.parm7")

os.chdir(OUTDIR)

print("=" * 80)
print("GB 방법으로 MMPBSA 재계산")
print("=" * 80)
print()
print(f"작업 디렉토리: {OUTDIR}")
print()

# Step 1: 물/이온 제거 좌표 확인 (이미 있으면 재사용)
if not (OUTDIR / "complex.rst7").exists():
    print("Step 1: 물/이온 제거 좌표 생성")
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
    
    if result.returncode != 0:
        print("❌ cpptraj 실패")
        print(result.stderr)
        exit(1)
    
    print("✅ 좌표 파일 생성: complex.rst7")
    print()
else:
    print("✅ 기존 좌표 파일 사용: complex.rst7")
    print()

# Step 2: GB MMPBSA 입력 파일 생성
print("Step 2: GB MMPBSA 입력 파일 생성")
print("-" * 80)

mmpbsa_gb_input = """MMPBSA calculation with GB (igb=8)
&general
startframe=1
endframe=1
interval=1
verbose=2
/

&gb
igb=8,
saltcon=0.150,
/
"""

with open("mmpbsa_gb.in", "w") as f:
    f.write(mmpbsa_gb_input)

print("✅ mmpbsa_gb.in 생성 완료")
print()

# Step 3: GB MMPBSA 실행
print("Step 3: GB MMPBSA 실행")
print("=" * 80)

cmd = [
    "MMPBSA.py",
    "-O",
    "-i", "mmpbsa_gb.in",
    "-o", "RESULTS_GB.dat",
    "-sp", "complex.prmtop",
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
    print("✅ GB MMPBSA 성공!")
    print("=" * 80)
    print()
    
    if Path("RESULTS_GB.dat").exists():
        print("GB 결과:")
        print("-" * 80)
        with open("RESULTS_GB.dat", "r") as f:
            content = f.read()
            print(content)
        
        # 핵심 결과 추출
        print()
        print("=" * 80)
        print("핵심 결과 비교")
        print("=" * 80)
        print()
        
        # GB 결과에서 DELTA 값 찾기
        for line in content.split('\n'):
            if 'DELTA G gas' in line or 'DELTA G solv' in line or 'DELTA TOTAL' in line:
                print(line)
        
        # PB 결과와 비교
        pb_file = Path("FINAL_RESULTS_MMPBSA.dat")
        if pb_file.exists():
            print()
            print("PB 결과 (참고):")
            with open(pb_file, "r") as f:
                pb_content = f.read()
                for line in pb_content.split('\n'):
                    if 'DELTA G gas' in line or 'DELTA G solv' in line or 'DELTA TOTAL' in line:
                        print(line)
    else:
        print("⚠️  결과 파일이 생성되지 않았습니다.")
else:
    print(f"❌ GB MMPBSA 실패 (exit code: {result.returncode})")
    print("=" * 80)
    exit(1)
