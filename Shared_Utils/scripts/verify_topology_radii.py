#!/usr/bin/env python3
"""
생성된 topology 파일들의 radii 확인
"""

import parmed as pmd
import sys

print("=" * 80)
print("Topology 파일 Radii 확인")
print("=" * 80)
print()

# 출력 디렉토리 읽기
with open('/tmp/mmpbsa_outdir.txt', 'r') as f:
    outdir = f.read().strip()

print(f"출력 디렉토리: {outdir}")
print()

topology_files = {
    'complex': f'{outdir}/complex.prmtop',
    'receptor': f'{outdir}/receptor.prmtop',
    'ligand': f'{outdir}/ligand.prmtop'
}

all_ok = True

for name, filepath in topology_files.items():
    print(f"📁 {name}.prmtop")
    print(f"   파일: {filepath}")
    
    try:
        parm = pmd.load_file(filepath)
        radii = [a.solvent_radius for a in parm.atoms]
        unique_radii = sorted(set(radii))
        none_count = sum(r is None or r == 0 for r in radii)
        
        print(f"   총 원자: {len(parm.atoms)}")
        print(f"   고유 radii 값: {len(unique_radii)}")
        print(f"   샘플 radii: {unique_radii[:5]}")
        print(f"   Radii 없는 원자: {none_count}")
        
        if none_count == 0 and len(unique_radii) > 1:
            print(f"   ✅ Radii 정상")
        else:
            print(f"   ❌ Radii 문제 발견")
            all_ok = False
    except Exception as e:
        print(f"   ❌ 에러: {e}")
        all_ok = False
    
    print()

print("=" * 80)
if all_ok:
    print("✅ 모든 topology 파일에 radii가 정상적으로 유지됨")
    print("   MMPBSA 계산 준비 완료!")
else:
    print("❌ 일부 topology 파일에 radii 문제 있음")
    sys.exit(1)
print("=" * 80)
