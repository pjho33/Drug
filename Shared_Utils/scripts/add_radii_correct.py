#!/usr/bin/env python3
"""
ParmEd 정석 방법으로 radii 추가
"""

import parmed as pmd
from parmed.tools.actions import changeRadii

print("=" * 80)
print("Adding mbondi2 radii to topology")
print("=" * 80)
print()

print("Loading original parm7...")
parm = pmd.load_file(
    '/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/GLUT1SDGComplex260110/amber/step5_input.parm7'
)
print(f"  Atoms: {len(parm.atoms)}")
print(f"  Residues: {len(parm.residues)}")
print()

print("Applying mbondi2 radii...")
action = changeRadii(parm, 'mbondi2')
action.execute()
print("  ✅ Radii applied")
print()

output_file = '/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/mmpbsa_amber/step5_input_radii.parm7'
print(f"Saving to: {output_file}")
parm.save(output_file, overwrite=True)
print("  ✅ Saved")
print()

# 파일 확인
import os
if os.path.exists(output_file):
    size_mb = os.path.getsize(output_file) / (1024 * 1024)
    print(f"✅ File created: {size_mb:.1f} MB")
else:
    print("❌ File not created!")

print()
print("=" * 80)
print("✅ Complete!")
print("=" * 80)
