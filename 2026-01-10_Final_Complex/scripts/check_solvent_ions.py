#!/usr/bin/env python3
"""
원본 topology에서 물/이온 residue 이름 확인
"""

import parmed as pmd
from collections import Counter

RADII_PARM = "/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/mmpbsa_amber/step5_input_radii.parm7"

print("=" * 80)
print("물/이온 Residue 확인")
print("=" * 80)
print()

parm = pmd.load_file(RADII_PARM)
res_types = Counter([r.name for r in parm.residues])

# 단백질/리간드 제외
protein_res = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
               'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
               'HIE', 'HID', 'HIP', 'CYX', 'HSD']
ligand_res = ['SDG']
lipid_res = ['POPC', 'POPE', 'POPS', 'PSM', 'CHL1']

water_ion_res = {}
for res_name, count in res_types.items():
    if res_name not in protein_res + ligand_res + lipid_res:
        water_ion_res[res_name] = count

print("물/이온 Residues:")
for res_name, count in sorted(water_ion_res.items()):
    print(f"  {res_name}: {count}")

print()
print("Strip mask 생성:")
strip_mask = ":" + ",".join(water_ion_res.keys())
print(f"  {strip_mask}")
print()

# 원자 수 계산
water_ion_atoms = sum(len(r.atoms) for r in parm.residues if r.name in water_ion_res)
total_atoms = len(parm.atoms)
complex_atoms = total_atoms - water_ion_atoms

print(f"총 원자: {total_atoms}")
print(f"물/이온 원자: {water_ion_atoms}")
print(f"Complex 예상 원자 (물/이온 제거 후): {complex_atoms}")
print()

print("=" * 80)
