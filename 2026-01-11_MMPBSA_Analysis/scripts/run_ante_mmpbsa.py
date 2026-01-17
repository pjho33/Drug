#!/usr/bin/env python3
"""
ante-MMPBSA.py 실행 (radii 포함 topology)
"""

import subprocess
import os
from pathlib import Path

MMPBSA_DIR = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/mmpbsa_amber")

# 현재 작업 디렉토리 확인 (호출자가 지정한 디렉토리 사용)
print("=" * 80)
print("ante-MMPBSA.py 실행 (radii 포함 topology)")
print("=" * 80)
print()
print(f"PWD: {os.getcwd()}")
print(f"SCRIPT: {__file__}")
print()

print("입력 파일:")
print("  step5_input_radii.parm7 (radii 포함)")
print()

print("출력 파일:")
print("  complex.prmtop")
print("  receptor.prmtop")
print("  ligand.prmtop")
print()

print("Masks:")
print("  Receptor: :1-305")
print("  Ligand: :306")
print("  Strip (물/이온): :POT,CLA,WAT")
print()

# 기존 출력 파일 삭제 (충돌 방지)
for fname in ["complex.prmtop", "receptor.prmtop", "ligand.prmtop"]:
    if os.path.exists(fname):
        os.remove(fname)
        print(f"  기존 파일 삭제: {fname}")

print()
print("=" * 80)
print("실행 중...")
print("=" * 80)
print()

cmd = [
    "ante-MMPBSA.py",
    "-p", str(MMPBSA_DIR / "step5_input_radii.parm7"),
    "-c", "complex.prmtop",
    "-r", "receptor.prmtop",
    "-l", "ligand.prmtop",
    "-s", ":POT,CLA,WAT",  # 물/이온만 strip
    "-n", ":1-305"  # receptor (ligand는 자동 = receptor 제외한 나머지)
]

result = subprocess.run(cmd, capture_output=True, text=True)

print(result.stdout)
if result.stderr:
    print(result.stderr)

print()
print("=" * 80)
if result.returncode == 0:
    print("✅ ante-MMPBSA.py 성공!")
    print("=" * 80)
    print()
    
    print("생성된 파일:")
    for fname in ["complex.prmtop", "receptor.prmtop", "ligand.prmtop"]:
        if Path(fname).exists():
            size_mb = Path(fname).stat().st_size / (1024 * 1024)
            print(f"  {fname}: {size_mb:.1f} MB")
    print()
    
    # 원자수 확인
    print("원자수 확인:")
    try:
        import parmed as pmd
        for fname in ["complex.prmtop", "receptor.prmtop", "ligand.prmtop"]:
            if Path(fname).exists():
                p = pmd.load_file(fname)
                print(f"  {fname}: {len(p.atoms)} atoms, {len(p.residues)} residues")
    except Exception as e:
        print(f"  원자수 확인 실패: {e}")
    print()
    
    print("다음 단계: MMPBSA.py 실행")
else:
    print(f"❌ ante-MMPBSA.py 실패 (exit code: {result.returncode})")
    print("=" * 80)
    exit(1)
