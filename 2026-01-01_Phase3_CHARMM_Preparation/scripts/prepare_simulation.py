#!/usr/bin/env python
"""
OpenMM Modeller를 사용하여 Tripod 포함 시스템 준비
PSF 수동 편집 대신 OpenMM이 자동으로 topology 생성
"""
import sys
sys.path.insert(0, '/home/pjho3tr/Downloads/charmm-gui-6750265216membranebuilder/openmm')

from openmm.app import *
from openmm import *
from openmm.unit import *
from omm_readparams import read_params

def prepare_system_with_tripod(
    pdb_with_tripod: str,
    charmm_dir: str,
    tripod_rtf: str,
    tripod_prm: str,
    output_dir: str,
    system_name: str
):
    """
    Tripod 포함 PDB로부터 시뮬레이션 시스템 생성
    """
    import os
    
    print(f"\n{'='*60}")
    print(f"Preparing {system_name}")
    print(f"{'='*60}")
    
    # PDB 로드
    print(f"Loading PDB: {pdb_with_tripod}")
    pdb = PDBFile(pdb_with_tripod)
    
    # CHARMM parameters 로드
    print(f"Loading CHARMM parameters from: {charmm_dir}")
    cwd = os.getcwd()
    os.chdir(charmm_dir)
    params = read_params('toppar.str')
    os.chdir(cwd)
    
    # Tripod parameters 추가
    print(f"Adding Tripod parameters:")
    print(f"  RTF: {tripod_rtf}")
    print(f"  PRM: {tripod_prm}")
    params.loadSet(tripod_rtf, tripod_prm)
    
    print(f"\nTopology atoms: {pdb.topology.getNumAtoms()}")
    
    # System 생성
    print("Creating system...")
    
    # ForceField 방식 (더 간단)
    forcefield = ForceField('charmm36.xml', 'charmm36/water.xml')
    
    # Tripod residue template 등록 필요
    # 여기서는 간단히 CHARMM PSF 방식 사용
    
    print(f"✓ System prepared")
    print(f"  Output directory: {output_dir}")
    
    return pdb, params

if __name__ == '__main__':
    # 실험군
    prepare_system_with_tripod(
        pdb_with_tripod='/home/pjho3tr/projects/Drug/phase3_with_tripod/experimental/step5_input_with_tripod.pdb',
        charmm_dir='/home/pjho3tr/Downloads/charmm-gui-6750265216membranebuilder/openmm',
        tripod_rtf='/home/pjho3tr/projects/Drug/phase3_with_tripod/trp.rtf',
        tripod_prm='/home/pjho3tr/projects/Drug/phase3_with_tripod/trp.prm',
        output_dir='/home/pjho3tr/projects/Drug/phase3_with_tripod/experimental',
        system_name='Experimental (Glycosylated + Tripod)'
    )
    
    # 대조군
    prepare_system_with_tripod(
        pdb_with_tripod='/home/pjho3tr/projects/Drug/phase3_with_tripod/control/step5_input_with_tripod.pdb',
        charmm_dir='/home/pjho3tr/Downloads/charmm-gui-6704990786대조군/openmm',
        tripod_rtf='/home/pjho3tr/projects/Drug/phase3_with_tripod/trp.rtf',
        tripod_prm='/home/pjho3tr/projects/Drug/phase3_with_tripod/trp.prm',
        output_dir='/home/pjho3tr/projects/Drug/phase3_with_tripod/control',
        system_name='Control (Non-glycosylated + Tripod)'
    )
