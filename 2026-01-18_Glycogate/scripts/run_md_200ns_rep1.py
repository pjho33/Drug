#!/usr/bin/env python3
"""
TRIS-PEG24-Lglucose 200ns MD 시뮬레이션 (Replica 1)
배경 실행용 스크립트
"""

import os
import sys
from pathlib import Path

def check_openmm():
    """OpenMM 설치 및 GPU 확인"""
    print("=" * 80)
    print("OpenMM 환경 확인")
    print("=" * 80)
    print()
    
    try:
        import openmm
        from openmm import app, unit
        print(f"✅ OpenMM version: {openmm.version.version}")
        print()
        
        # Platform 확인
        print("사용 가능한 Platform:")
        for i in range(openmm.Platform.getNumPlatforms()):
            platform = openmm.Platform.getPlatform(i)
            print(f"  - {platform.getName()}")
        print()
        
        # CUDA 확인
        try:
            cuda = openmm.Platform.getPlatformByName('CUDA')
            print(f"✅ CUDA Platform 사용 가능")
            print()
        except Exception:
            print("⚠️  CUDA Platform 없음 (CPU 사용)")
            print()
        
        return True
        
    except ImportError:
        print("❌ OpenMM이 설치되지 않았습니다.")
        return False


def setup_simulation(pdb_file, psf_file, toppar_dir, output_prefix, 
                     temperature=300, steps=100000000, replica=1):
    """OpenMM 시뮬레이션 설정"""
    
    from openmm import app, unit, Platform
    from openmm.app import PDBFile, CharmmPsfFile, CharmmParameterSet, CharmmCrdFile
    import openmm
    
    print("=" * 80)
    print(f"OpenMM 시뮬레이션 설정 (Replica {replica})")
    print("=" * 80)
    print()
    
    # 1. 구조 및 Topology 로드
    print("Step 1: 구조 및 Topology 로드")
    print("-" * 80)
    
    psf = CharmmPsfFile(psf_file)
    pdb = PDBFile(pdb_file)
    
    # CRD 파일에서 box 정보 읽기
    crd_file = pdb_file.replace('.pdb', '.crd')
    if os.path.exists(crd_file):
        crd = CharmmCrdFile(crd_file)
        if crd.boxVectors is not None:
            print(f"✅ CRD에서 box vectors 로드")
            psf.setBox(crd.boxVectors[0][0], crd.boxVectors[1][1], crd.boxVectors[2][2])
    
    # 또는 crystal_image.str에서 box 정보 읽기
    crystal_file = str(Path(toppar_dir).parent / "crystal_image.str")
    if os.path.exists(crystal_file):
        with open(crystal_file, 'r') as f:
            for line in f:
                if 'SET A' in line or 'SET B' in line or 'SET C' in line:
                    print(f"Box info: {line.strip()}")
                    # A, B, C 값 추출하여 box 설정
                    if 'SET A' in line:
                        a_val = float(line.split('=')[1].strip())
                    elif 'SET B' in line:
                        b_val = float(line.split('=')[1].strip())
                    elif 'SET C' in line:
                        c_val = float(line.split('=')[1].strip())
        try:
            psf.setBox(a_val*unit.angstrom, b_val*unit.angstrom, c_val*unit.angstrom)
            print(f"✅ crystal_image.str에서 box 설정: {a_val} x {b_val} x {c_val} Å")
        except:
            pass
    
    print(f"✅ PSF: {psf_file}")
    print(f"✅ PDB: {pdb_file}")
    print(f"   원자 수: {psf.topology.getNumAtoms()}")
    print()
    
    # 2. Force field 로드
    print("Step 2: Force field 로드")
    print("-" * 80)
    
    # CHARMM36 파라미터 - 순서대로 로드 (RTF -> PRM -> STR)
    import glob
    from pathlib import Path
    
    # 1. RTF 파일 (topology) 먼저
    rtf_files = sorted(glob.glob(f"{toppar_dir}/*.rtf"))
    # 2. PRM 파일 (parameters)
    prm_files = sorted(glob.glob(f"{toppar_dir}/*.prm"))
    # 3. STR 파일 (stream files)
    str_files = sorted(glob.glob(f"{toppar_dir}/*.str"))
    
    # 4. 리간드 파라미터 추가 (Ligand Reader/lig)
    ligand_dir = Path(toppar_dir).parent.parent / "Ligand reader" / "lig"
    ligand_rtf = str(ligand_dir / "lig.rtf")
    ligand_prm = str(ligand_dir / "lig.prm")
    
    if Path(ligand_rtf).exists():
        rtf_files.append(ligand_rtf)
        print(f"✅ 리간드 RTF 추가: {ligand_rtf}")
    if Path(ligand_prm).exists():
        prm_files.append(ligand_prm)
        print(f"✅ 리간드 PRM 추가: {ligand_prm}")
    
    toppar_files = rtf_files + prm_files + str_files
    
    print(f"Loading {len(toppar_files)} parameter files...")
    print(f"  RTF: {len(rtf_files)}, PRM: {len(prm_files)}, STR: {len(str_files)}")
    
    params = CharmmParameterSet(*toppar_files)
    
    print("✅ CHARMM36 force field 로드 완료")
    print()
    
    # 3. System 생성
    print("Step 3: System 생성")
    print("-" * 80)
    
    system = psf.createSystem(
        params,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.2*unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
        ewaldErrorTolerance=0.0005
    )
    
    print("✅ System 생성 완료")
    print()
    
    # 4. Integrator 설정
    print("Step 4: Integrator 설정")
    print("-" * 80)
    
    integrator = openmm.LangevinMiddleIntegrator(
        temperature*unit.kelvin,
        1.0/unit.picosecond,
        0.002*unit.picoseconds
    )
    
    print(f"✅ Langevin integrator")
    print(f"   Temperature: {temperature} K")
    print(f"   Friction: 1.0 ps^-1")
    print(f"   Timestep: 2 fs")
    print()
    
    # 5. Barostat 추가 (NPT)
    print("Step 5: Barostat 추가 (NPT)")
    print("-" * 80)
    
    system.addForce(openmm.MonteCarloBarostat(
        1.0*unit.bar,
        temperature*unit.kelvin,
        25
    ))
    
    print("✅ Monte Carlo barostat")
    print("   Pressure: 1.0 bar")
    print()
    
    # 6. Platform 선택
    print("Step 6: Platform 선택")
    print("-" * 80)
    
    try:
        platform = Platform.getPlatformByName('CUDA')
        properties = {'Precision': 'mixed'}
        print("✅ CUDA Platform (mixed precision)")
    except Exception:
        try:
            platform = Platform.getPlatformByName('OpenCL')
            properties = {}
            print("✅ OpenCL Platform")
        except Exception:
            platform = Platform.getPlatformByName('CPU')
            properties = {}
            print("⚠️  CPU Platform (느림)")
    print()
    
    # 7. Simulation 생성
    print("Step 7: Simulation 생성")
    print("-" * 80)
    
    simulation = app.Simulation(
        psf.topology,
        system,
        integrator,
        platform,
        properties
    )
    
    simulation.context.setPositions(pdb.positions)
    
    print("✅ Simulation 생성 완료")
    print()
    
    # 8. Reporters 설정
    print("Step 8: Reporters 설정")
    print("-" * 80)
    
    # DCD trajectory (10 ps 간격)
    simulation.reporters.append(
        app.DCDReporter(f"{output_prefix}.dcd", 5000)
    )
    
    # State data (에너지, 온도 등)
    simulation.reporters.append(
        app.StateDataReporter(
            f"{output_prefix}.log",
            5000,
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            volume=True,
            density=True,
            speed=True
        )
    )
    
    # Checkpoint (100 ps 간격)
    simulation.reporters.append(
        app.CheckpointReporter(f"{output_prefix}.chk", 50000)
    )
    
    print(f"✅ DCD: {output_prefix}.dcd (10 ps 간격)")
    print(f"✅ LOG: {output_prefix}.log")
    print(f"✅ CHK: {output_prefix}.chk (100 ps 간격)")
    print()
    
    return simulation


def run_equilibration(simulation, output_prefix):
    """에너지 최소화 및 평형화"""
    
    from openmm import unit
    
    print("=" * 80)
    print("에너지 최소화 및 평형화")
    print("=" * 80)
    print()
    
    # 1. 에너지 최소화
    print("Step 1: 에너지 최소화")
    print("-" * 80)
    
    print("최소화 중...")
    simulation.minimizeEnergy(maxIterations=1000)
    
    state = simulation.context.getState(getEnergy=True)
    print(f"✅ 최소화 완료")
    print(f"   Potential Energy: {state.getPotentialEnergy()}")
    print()
    
    # 2. NVT 평형화 (100 ps)
    print("Step 2: NVT 평형화 (100 ps)")
    print("-" * 80)
    
    simulation.context.setVelocitiesToTemperature(300)
    print("NVT 평형화 중...")
    simulation.step(50000)  # 100 ps
    
    print("✅ NVT 평형화 완료")
    print()
    
    # 3. 초기 구조 저장
    print("Step 3: 평형화 후 구조 저장")
    print("-" * 80)
    
    from openmm.app import PDBFile
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(
        simulation.topology,
        positions,
        open(f"{output_prefix}_equilibrated.pdb", 'w')
    )
    
    print(f"✅ {output_prefix}_equilibrated.pdb")
    print()


def run_production(simulation, steps, output_prefix):
    """Production MD 실행"""
    
    print("=" * 80)
    print("Production MD 실행")
    print("=" * 80)
    print()
    
    print(f"총 스텝: {steps:,}")
    print(f"시뮬레이션 시간: {steps * 0.002 / 1000:.1f} ns")
    print()
    
    print("실행 중...")
    print("진행 상황은 로그 파일을 확인하세요.")
    print()
    
    simulation.step(steps)
    
    print()
    print("✅ Production MD 완료!")
    print()


def main():
    """메인 함수"""
    
    # OpenMM 확인
    if not check_openmm():
        sys.exit(1)
    
    # 경로 설정
    work_dir = Path("/home/pjho3tr/projects/Drug/2026-01-18_Glycogate")
    charmm_gui_dir = work_dir / "data" / "solution builder" / "openmm"
    results_dir = work_dir / "results" / "md_200ns_rep1"
    
    results_dir.mkdir(parents=True, exist_ok=True)
    
    # 입력 파일
    pdb_file = str(charmm_gui_dir / "step3_input.pdb")
    psf_file = str(charmm_gui_dir / "step3_input.psf")
    toppar_dir = str(work_dir / "data" / "solution builder" / "toppar")
    
    # 출력 prefix
    output_prefix = str(results_dir / "md_rep1")
    
    # 시뮬레이션 설정
    simulation = setup_simulation(
        pdb_file, psf_file, toppar_dir, output_prefix,
        temperature=300,
        steps=100000000,  # 200 ns
        replica=1
    )
    
    # 평형화
    run_equilibration(simulation, output_prefix)
    
    # Production MD
    run_production(simulation, 100000000, output_prefix)  # 200 ns
    
    print("=" * 80)
    print("✅ 모든 작업 완료!")
    print("=" * 80)
    print()
    print(f"출력 파일: {results_dir}/md_rep1.*")
    print()


if __name__ == "__main__":
    main()
