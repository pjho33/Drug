#!/usr/bin/env python3
"""
원본 parm7 파일에 mbondi2 radii 추가
"""

import parmed as pmd
from pathlib import Path

# 경로 설정
AMBER_DIR = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/GLUT1SDGComplex260110/amber")
MMPBSA_DIR = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/mmpbsa_amber")

original_parm = AMBER_DIR / "step5_input.parm7"
output_parm = MMPBSA_DIR / "step5_input_radii.parm7"

print("=" * 80)
print("Adding Radii to Topology File")
print("=" * 80)
print()

# 원본 파일 로드
print(f"Loading: {original_parm}")
parm = pmd.load_file(str(original_parm))
print(f"  Atoms: {len(parm.atoms)}")
print(f"  Residues: {len(parm.residues)}")
print()

# radii 추가
print("Adding mbondi2 radii...")
parm = pmd.tools.changeRadii(parm, 'mbondi2')
print("  ✅ Radii added")
print()

# 저장
print(f"Saving to: {output_parm}")
parm.save(str(output_parm), overwrite=True)
print("  ✅ Saved")
print()

# 확인
if output_parm.exists():
    size_mb = output_parm.stat().st_size / (1024 * 1024)
    print(f"✅ File created: {output_parm.name} ({size_mb:.1f} MB)")
else:
    print("❌ File not created!")
    exit(1)

print()
print("=" * 80)
print("✅ Complete!")
print("=" * 80)
