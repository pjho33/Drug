#!/usr/bin/env python3
"""
complex.prmtop에 ligand가 포함되어 있는지 확인
"""

import parmed as pmd
from pathlib import Path

# 출력 디렉토리 읽기
with open('/tmp/mmpbsa_outdir.txt', 'r') as f:
    outdir = Path(f.read().strip())

complex_prmtop = outdir / "complex.prmtop"

print("=" * 80)
print("Complex Topology 확인")
print("=" * 80)
print()
print(f"파일: {complex_prmtop}")
print()

parm = pmd.load_file(str(complex_prmtop))
resnames = sorted(set(r.name for r in parm.residues))

print(f"총 residues: {len(parm.residues)}")
print(f"총 atoms: {len(parm.atoms)}")
print()

has_sdg = "SDG" in resnames
print(f"Has SDG? {has_sdg}")
print()

if has_sdg:
    sdg_residues = [r for r in parm.residues if r.name == "SDG"]
    print(f"✅ SDG 발견! ({len(sdg_residues)}개)")
    for r in sdg_residues:
        print(f"   Residue {r.idx+1}: {r.name}, {len(r.atoms)} atoms")
else:
    print("❌ SDG 없음 - ligand가 제거되었습니다!")

print()
print("마지막 20개 residue names:")
print([r.name for r in parm.residues[-20:]])
print()

print("=" * 80)
if has_sdg:
    print("✅ Complex에 ligand 포함됨 - 정상")
else:
    print("❌ Complex에서 ligand 제거됨 - ante-MMPBSA mask 수정 필요")
print("=" * 80)
