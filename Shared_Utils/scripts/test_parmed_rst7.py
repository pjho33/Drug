#!/usr/bin/env python3
"""
ParmEd rst7 인식 테스트
"""

from parmed import load_file
import sys

parm_path = "/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/mmpbsa_amber/step5_input_radii.parm7"
rst_path = "/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/GLUT1SDGComplex260110/amber/step5_input.rst7"

print("=" * 80)
print("ParmEd rst7 인식 테스트")
print("=" * 80)
print()

print("Loading parm...")
try:
    parm = load_file(parm_path)
    parm_natom = parm.ptr('natom')
    print(f"✅ parm natom: {parm_natom}")
except Exception as e:
    print(f"❌ parm 로딩 실패: {e}")
    sys.exit(1)

print()
print("Loading rst7...")
try:
    rst = load_file(rst_path)
    rst_natom = rst.coordinates.shape[0]
    print(f"✅ rst7 natom: {rst_natom}")
    print()
    
    if parm_natom == rst_natom:
        print("✅ 원자 수 일치 - Minimization 진행 가능")
    else:
        print(f"❌ 원자 수 불일치 (차이: {abs(parm_natom - rst_natom)})")
        
except Exception as e:
    print(f"❌ rst7 로딩 실패: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print()
print("=" * 80)
