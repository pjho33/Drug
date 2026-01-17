#!/usr/bin/env python
"""
PDB 포맷 수정: 원자 번호가 99999를 초과하면 정렬이 깨지므로 재번호화
"""
import sys

def fix_pdb_format(input_pdb, output_pdb):
    """PDB 원자 번호를 1부터 재번호화"""
    with open(input_pdb, 'r') as f:
        lines = f.readlines()
    
    atom_num = 1
    output_lines = []
    
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            # 원자 번호 재할당 (1-5 컬럼)
            record = line[:6]
            rest = line[11:]  # 원자 이름부터 끝까지
            new_line = f"{record}{atom_num:5d}{rest}"
            output_lines.append(new_line)
            atom_num += 1
        else:
            output_lines.append(line)
    
    with open(output_pdb, 'w') as f:
        f.writelines(output_lines)
    
    print(f"Fixed: {input_pdb} -> {output_pdb}")
    print(f"  Total atoms: {atom_num - 1}")

# 실험군
fix_pdb_format(
    '/home/pjho3tr/projects/Drug/phase3_with_tripod/experimental/step5_input_with_tripod.pdb',
    '/home/pjho3tr/projects/Drug/phase3_with_tripod/experimental/step5_input_with_tripod_fixed.pdb'
)

# 대조군
fix_pdb_format(
    '/home/pjho3tr/projects/Drug/phase3_with_tripod/control/step5_input_with_tripod.pdb',
    '/home/pjho3tr/projects/Drug/phase3_with_tripod/control/step5_input_with_tripod_fixed.pdb'
)

print("\n✓ PDB files fixed!")
