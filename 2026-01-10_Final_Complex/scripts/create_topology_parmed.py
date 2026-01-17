#!/usr/bin/env python3
"""
ParmEd로 SDG 기준 topology 분리
"""

import parmed as pmd
from pathlib import Path
import os

RADII_PARM = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/mmpbsa_amber/step5_input_radii.parm7")
OUTDIR = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/mmpbsa_amber/parmed_out")

OUTDIR.mkdir(exist_ok=True)
os.chdir(OUTDIR)

print("=" * 80)
print("ParmEd로 Topology 분리 (SDG 기준)")
print("=" * 80)
print()

# 원본 로드
print("원본 topology 로드...")
parm = pmd.load_file(str(RADII_PARM))
print(f"  총 atoms: {len(parm.atoms)}")
print(f"  총 residues: {len(parm.residues)}")
print()

# SDG residue 확인
sdg_residues = [r for r in parm.residues if r.name == "SDG"]
print(f"SDG residues: {len(sdg_residues)}")
if sdg_residues:
    for r in sdg_residues:
        print(f"  Residue {r.idx+1}: {r.name}, {len(r.atoms)} atoms")
print()

# Step 1: Complex (물/이온 제거)
print("Step 1: Complex topology 생성 (물/이온 제거)")
print("-" * 80)

complex_parm = pmd.load_file(str(RADII_PARM))
# 물/이온 제거
water_ion_mask = pmd.amber.AmberMask(complex_parm, ":POT,CLA,WAT")
complex_parm.strip(water_ion_mask)

complex_file = "complex.prmtop"
complex_parm.save(complex_file, overwrite=True)
print(f"✅ {complex_file}: {len(complex_parm.atoms)} atoms, {len(complex_parm.residues)} residues")
print()

# Step 2: Receptor (물/이온 + SDG 제거)
print("Step 2: Receptor topology 생성 (물/이온 + SDG 제거)")
print("-" * 80)

receptor_parm = pmd.load_file(str(RADII_PARM))
# 물/이온 + SDG 제거
strip_mask = pmd.amber.AmberMask(receptor_parm, ":POT,CLA,WAT,SDG")
receptor_parm.strip(strip_mask)

receptor_file = "receptor.prmtop"
receptor_parm.save(receptor_file, overwrite=True)
print(f"✅ {receptor_file}: {len(receptor_parm.atoms)} atoms, {len(receptor_parm.residues)} residues")
print()

# Step 3: Ligand (SDG만 남김)
print("Step 3: Ligand topology 생성 (SDG만)")
print("-" * 80)

ligand_parm = pmd.load_file(str(RADII_PARM))
# SDG 제외한 모든 것 제거
keep_sdg_mask = pmd.amber.AmberMask(ligand_parm, "!:SDG")
ligand_parm.strip(keep_sdg_mask)

ligand_file = "ligand.prmtop"
ligand_parm.save(ligand_file, overwrite=True)
print(f"✅ {ligand_file}: {len(ligand_parm.atoms)} atoms, {len(ligand_parm.residues)} residues")
print()

# 검증
print("=" * 80)
print("검증")
print("=" * 80)

complex_atoms = len(complex_parm.atoms)
receptor_atoms = len(receptor_parm.atoms)
ligand_atoms = len(ligand_parm.atoms)

print(f"Complex atoms: {complex_atoms}")
print(f"Receptor atoms: {receptor_atoms}")
print(f"Ligand atoms: {ligand_atoms}")
print(f"Receptor + Ligand: {receptor_atoms + ligand_atoms}")
print()

# Ligand에 SDG만 있는지 확인
ligand_resnames = set(r.name for r in ligand_parm.residues)
print(f"Ligand residue names: {ligand_resnames}")

if ligand_resnames == {"SDG"}:
    print("✅ Ligand에 SDG만 포함됨 - 정상")
else:
    print("⚠️  Ligand에 SDG 외 residue 포함됨")

print()
print("=" * 80)
print(f"✅ Topology 생성 완료: {OUTDIR}")
print("=" * 80)

# 출력 디렉토리 저장
with open('/tmp/mmpbsa_outdir.txt', 'w') as f:
    f.write(str(OUTDIR))
