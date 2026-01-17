#!/usr/bin/env python3
"""
Radii 적용 확인
"""

import parmed as pmd

print("=" * 80)
print("Verifying radii in topology file")
print("=" * 80)
print()

parm = pmd.load_file(
    '/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/mmpbsa_amber/step5_input_radii.parm7'
)

radii = [a.solvent_radius for a in parm.atoms]
unique_radii = sorted(set(radii))
none_count = sum(r is None or r == 0 for r in radii)

print(f"Total atoms: {len(parm.atoms)}")
print(f"Unique radii values: {len(unique_radii)}")
print(f"Sample radii: {unique_radii[:10]}")
print(f"Atoms without radii: {none_count}")
print()

if none_count == 0 and len(unique_radii) > 1:
    print("✅ Radii successfully applied!")
    print("   Ready for MMPBSA calculation")
else:
    print("❌ Radii not properly applied")
    exit(1)

print()
print("=" * 80)
