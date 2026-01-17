#!/usr/bin/env python3
"""
ligand.prmtop의 residue 구성 확인
"""

import parmed as pmd
from collections import Counter
from pathlib import Path

# 출력 디렉토리 읽기
with open('/tmp/mmpbsa_outdir.txt', 'r') as f:
    outdir = Path(f.read().strip())

ligand_prmtop = outdir / "ligand.prmtop"

print("=" * 80)
print("Ligand Topology Residue 구성 확인")
print("=" * 80)
print()
print(f"파일: {ligand_prmtop}")
print()

p = pmd.load_file(str(ligand_prmtop))
c = Counter(r.name for r in p.residues)

print(f"총 atoms: {len(p.atoms)}")
print(f"총 residues: {len(p.residues)}")
print()

print("Residue 구성 (상위 20개):")
for name, n in c.most_common(20):
    print(f"  {name}: {n}")

print()
print("=" * 80)

# SDG 확인
has_sdg = "SDG" in c
lipids = ['POPC', 'POPE', 'POPS', 'PSM', 'CHL1']
has_lipids = any(lipid in c for lipid in lipids)

if has_sdg and not has_lipids:
    print("✅ Ligand에 SDG만 포함됨 - 정상")
elif has_sdg and has_lipids:
    print("❌ Ligand에 SDG + lipids 혼입됨 - 분리 방식 수정 필요")
    print(f"   혼입된 lipids: {[l for l in lipids if l in c]}")
else:
    print("❌ Ligand에 SDG 없음 - 문제 있음")

print("=" * 80)
