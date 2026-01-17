#!/usr/bin/env python3
"""
Topology와 좌표 파일의 원자 수 일치 확인
"""

from parmed import load_file

PARM = "/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/mmpbsa_amber/step5_input_radii.parm7"
RST7 = "/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/GLUT1SDGComplex260110/amber/step5_input.rst7"

print("=" * 80)
print("Topology/좌표 원자 수 확인")
print("=" * 80)
print()

parm = load_file(PARM)
rst = load_file(RST7)

parm_atoms = parm.ptr("natom")
rst_atoms = rst.coordinates.shape[0]

print(f"Topology atoms: {parm_atoms}")
print(f"RST7 atoms:     {rst_atoms}")
print()

if parm_atoms == rst_atoms:
    print("✅ 원자 수 일치 - Minimization 진행 가능")
else:
    print("❌ 원자 수 불일치 - 문제 있음")
    print(f"   차이: {abs(parm_atoms - rst_atoms)} atoms")

print("=" * 80)
