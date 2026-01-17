#!/usr/bin/env python3
"""
gmx_MMPBSA 실행 (CHARMM36 force field)
"""

import subprocess
import os
from pathlib import Path

GROMACS_DIR = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/GLUT1SDGComplex260110/gromacs")
OUTDIR = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/gmx_mmpbsa_results")

OUTDIR.mkdir(parents=True, exist_ok=True)
os.chdir(OUTDIR)

print("=" * 80)
print("gmx_MMPBSA 실행 (CHARMM36)")
print("=" * 80)
print()
print(f"GROMACS 디렉토리: {GROMACS_DIR}")
print(f"출력 디렉토리: {OUTDIR}")
print()

# Step 1: 새 index 파일 생성 (Protein과 Ligand 그룹)
print("Step 1: Index 파일 생성")
print("-" * 80)

# gmx make_ndx로 그룹 생성
make_ndx_input = """
# Protein (residue 1-305 + membrane lipids 307-806)
r 1-305 | r 307-806
name 20 Receptor

# Ligand (residue 306 - SDG)
r 306
name 21 Ligand

q
"""

with open("make_ndx.txt", "w") as f:
    f.write(make_ndx_input)

cmd = [
    "gmx", "make_ndx",
    "-f", str(GROMACS_DIR / "step5_input.gro"),
    "-o", "index_mmpbsa.ndx"
]

print("명령:", " ".join(cmd))
result = subprocess.run(
    cmd,
    input=make_ndx_input,
    capture_output=True,
    text=True
)

if result.returncode != 0:
    print("❌ make_ndx 실패")
    print(result.stdout)
    print(result.stderr)
    exit(1)

print("✅ index_mmpbsa.ndx 생성 완료")
print()

# Step 2: gmx_MMPBSA 입력 파일 생성
print("Step 2: gmx_MMPBSA 입력 파일 생성")
print("-" * 80)

mmpbsa_input = """# gmx_MMPBSA input file
# CHARMM36 force field

&general
startframe=1,
endframe=1,
interval=1,
forcefields="charmm36",
PBRadii=4,
/

&gb
igb=8,
saltcon=0.150,
/
"""

with open("mmpbsa.in", "w") as f:
    f.write(mmpbsa_input)

print("✅ mmpbsa.in 생성 완료")
print()

# Step 3: gmx_MMPBSA 실행
print("Step 3: gmx_MMPBSA 실행")
print("=" * 80)

# PDB 파일 사용 (gmx_MMPBSA는 TPR 또는 PDB 필요)
pdb_file = GROMACS_DIR / "step5_input.pdb"
if not pdb_file.exists():
    print("❌ PDB 파일 없음")
    exit(1)

top_file = GROMACS_DIR / "topol.top"
if not top_file.exists():
    print("❌ topology 파일 없음")
    exit(1)

cmd = [
    "gmx_MMPBSA",
    "-O",
    "-i", "mmpbsa.in",
    "-cs", str(pdb_file),
    "-cp", str(top_file),
    "-ci", "index_mmpbsa.ndx",
    "-cg", "20", "21",  # Receptor, Ligand
    "-o", "FINAL_RESULTS_MMPBSA.dat",
    "-eo", "FINAL_RESULTS_MMPBSA.csv"
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
    print("✅ gmx_MMPBSA 성공!")
    print("=" * 80)
    print()
    
    # 결과 파일 확인
    result_files = [
        "FINAL_RESULTS_MMPBSA.dat",
        "FINAL_RESULTS_MMPBSA.csv",
        "_GMXMMPBSA_info"
    ]
    
    print("생성된 파일:")
    for fname in result_files:
        fpath = Path(fname)
        if fpath.exists():
            print(f"  ✅ {fname}")
        else:
            print(f"  ❌ {fname} (없음)")
    
    print()
    
    # 결과 출력
    if Path("FINAL_RESULTS_MMPBSA.dat").exists():
        print("결과:")
        print("-" * 80)
        with open("FINAL_RESULTS_MMPBSA.dat", "r") as f:
            print(f.read())
else:
    print(f"❌ gmx_MMPBSA 실패 (exit code: {result.returncode})")
    print("=" * 80)
    exit(1)
