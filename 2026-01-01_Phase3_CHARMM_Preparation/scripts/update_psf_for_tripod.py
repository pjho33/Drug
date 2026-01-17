#!/usr/bin/env python
"""
기존 PSF에 Tripod topology 추가
"""
import sys

def update_psf_with_tripod(
    original_psf: str,
    tripod_rtf: str,
    output_psf: str,
    n_tripod_atoms: int = 48
):
    """
    PSF에 Tripod 원자 추가
    
    Note: 간단한 버전 - Tripod을 별도 segment로 추가
    """
    print(f"Reading original PSF: {original_psf}")
    with open(original_psf, 'r') as f:
        psf_lines = f.readlines()
    
    # PSF 섹션 찾기
    natom_idx = None
    nbond_idx = None
    
    for i, line in enumerate(psf_lines):
        if '!NATOM' in line:
            natom_idx = i
        elif '!NBOND' in line:
            nbond_idx = i
            break
    
    if natom_idx is None:
        raise ValueError("Could not find !NATOM section")
    
    # 원래 원자 수
    original_natom = int(psf_lines[natom_idx].split()[0])
    print(f"  Original atoms: {original_natom}")
    
    # Tripod 원자 추가 (간단 버전)
    # 실제로는 trp.rtf를 파싱해야 하지만, 여기서는 placeholder
    tripod_atoms = []
    for i in range(n_tripod_atoms):
        atom_id = original_natom + i + 1
        # Placeholder: 실제 atom type은 trp.rtf에서 가져와야 함
        tripod_atoms.append(
            f"{atom_id:8d} TRP  1    TRP  C{i+1:<4s} C      0.000000       12.0110           0\n"
        )
    
    # 새 NATOM 라인
    new_natom = original_natom + n_tripod_atoms
    new_natom_line = f"{new_natom:8d} !NATOM\n"
    
    # PSF 재구성
    output_lines = (
        psf_lines[:natom_idx] +
        [new_natom_line] +
        psf_lines[natom_idx+1:nbond_idx] +
        tripod_atoms +
        psf_lines[nbond_idx:]
    )
    
    print(f"Writing updated PSF: {output_psf}")
    with open(output_psf, 'w') as f:
        f.writelines(output_lines)
    
    print(f"  New total atoms: {new_natom}")
    print(f"  ✓ PSF updated")

if __name__ == '__main__':
    print("="*60)
    print("WARNING: PSF 수동 업데이트는 복잡합니다")
    print("대신 OpenMM Modeller를 사용하는 것을 권장합니다")
    print("="*60)
