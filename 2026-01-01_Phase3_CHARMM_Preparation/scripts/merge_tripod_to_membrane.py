#!/usr/bin/env python
"""
기존 CHARMM-GUI membrane 시스템에 Tripod를 추가
- 기존: step5_input.pdb (GLUT1 + membrane + water + ions)
- 추가: Tripod 좌표 (복합체에서 추출)
- 결과: GLUT1 + Tripod + membrane + water + ions
"""
import sys
from pathlib import Path

def merge_tripod_to_membrane(
    charmm_pdb: str,
    complex_pdb: str,
    output_pdb: str
):
    """
    기존 CHARMM-GUI PDB에 Tripod 추가
    
    Args:
        charmm_pdb: 기존 CHARMM-GUI step5_input.pdb (membrane 포함)
        complex_pdb: Tripod 포함 복합체 PDB
        output_pdb: 출력 PDB (GLUT1 + Tripod + membrane)
    """
    print(f"Reading CHARMM-GUI system: {charmm_pdb}")
    with open(charmm_pdb, 'r') as f:
        charmm_lines = f.readlines()
    
    print(f"Reading complex (for Tripod): {complex_pdb}")
    with open(complex_pdb, 'r') as f:
        complex_lines = f.readlines()
    
    # Tripod 좌표 추출 (HETATM, residue name TRP)
    tripod_lines = []
    for line in complex_lines:
        if line.startswith('HETATM') and 'TRP' in line:
            tripod_lines.append(line)
    
    print(f"  Found {len(tripod_lines)} Tripod atoms")
    
    # 기존 CHARMM PDB에서 단백질 마지막 원자 번호 찾기
    last_protein_atom = 0
    protein_end_idx = 0
    
    for i, line in enumerate(charmm_lines):
        if line.startswith('ATOM'):
            last_protein_atom = i
    
    protein_end_idx = last_protein_atom + 1
    
    print(f"  Protein ends at line {protein_end_idx}")
    
    # Tripod atom 번호 재할당
    if charmm_lines[last_protein_atom].startswith('ATOM'):
        last_atom_num = int(charmm_lines[last_protein_atom][6:11].strip())
    else:
        last_atom_num = 0
    
    print(f"  Last protein atom number: {last_atom_num}")
    
    # Tripod 원자 번호 업데이트
    renumbered_tripod = []
    for i, line in enumerate(tripod_lines):
        new_atom_num = last_atom_num + i + 1
        new_line = f"HETATM{new_atom_num:5d}" + line[11:]
        renumbered_tripod.append(new_line)
    
    # 병합: 단백질 + Tripod + 나머지 (membrane, water, ions)
    output_lines = (
        charmm_lines[:protein_end_idx] +
        renumbered_tripod +
        charmm_lines[protein_end_idx:]
    )
    
    # 저장
    print(f"Writing merged system: {output_pdb}")
    with open(output_pdb, 'w') as f:
        f.writelines(output_lines)
    
    # 통계
    n_atoms_original = sum(1 for l in charmm_lines if l.startswith(('ATOM', 'HETATM')))
    n_atoms_new = sum(1 for l in output_lines if l.startswith(('ATOM', 'HETATM')))
    
    print(f"\n✓ Merge complete!")
    print(f"  Original atoms: {n_atoms_original}")
    print(f"  Added Tripod atoms: {len(tripod_lines)}")
    print(f"  Total atoms: {n_atoms_new}")
    print(f"  Output: {output_pdb}")

if __name__ == '__main__':
    # 실험군 (Glycosylated)
    print("="*60)
    print("EXPERIMENTAL GROUP (Glycosylated + Tripod)")
    print("="*60)
    merge_tripod_to_membrane(
        charmm_pdb='/home/pjho3tr/Downloads/charmm-gui-6750265216membranebuilder/openmm/step5_input.pdb',
        complex_pdb='/home/pjho3tr/projects/Drug/scripts/glut1_tripod_complex.pdb',
        output_pdb='/home/pjho3tr/projects/Drug/phase3_with_tripod/experimental/step5_input_with_tripod.pdb'
    )
    
    print("\n" + "="*60)
    print("CONTROL GROUP (Non-glycosylated + Tripod)")
    print("="*60)
    merge_tripod_to_membrane(
        charmm_pdb='/home/pjho3tr/Downloads/charmm-gui-6704990786대조군/openmm/step5_input.pdb',
        complex_pdb='/home/pjho3tr/projects/Drug/scripts/glut1_tripod_complex_control.pdb',
        output_pdb='/home/pjho3tr/projects/Drug/phase3_with_tripod/control/step5_input_with_tripod.pdb'
    )
    
    print("\n" + "="*60)
    print("✓ Both systems ready for MD simulation!")
    print("="*60)
