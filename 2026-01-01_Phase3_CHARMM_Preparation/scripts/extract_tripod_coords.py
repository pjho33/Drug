#!/usr/bin/env python
"""
Phase 2 결과에서 Tripod 좌표를 추출하여 CHARMM-GUI에서 사용할 수 있는 형식으로 변환
"""
import MDAnalysis as mda

phase2_pdb = '/home/pjho3tr/projects/Drug/results/phase2_rep1/prod_tripod_rep1_final.pdb'

u = mda.Universe(phase2_pdb)
ligand = u.select_atoms('not protein and not resname HOH NA CL K TIP3 WAT SOL')

print("Tripod 좌표 정보:")
print(f"  원자 수: {len(ligand)}")
print(f"  중심 좌표: {ligand.center_of_mass()}")
print(f"\nBinding pocket 근처 단백질 잔기:")

protein = u.select_atoms('protein')
nearby = protein.select_atoms(f'around 5 (resname {ligand.residues[0].resname})')
print(f"  5 Å 이내 잔기: {len(nearby.residues)}개")
for res in nearby.residues[:10]:
    print(f"    {res.resname}{res.resid}")

# Tripod만 포함된 PDB 저장
ligand.write('tripod_only.pdb')
print(f"\n✓ Tripod 좌표 저장: tripod_only.pdb")
print(f"  이 파일을 CHARMM-GUI에서 ligand position reference로 사용 가능")
