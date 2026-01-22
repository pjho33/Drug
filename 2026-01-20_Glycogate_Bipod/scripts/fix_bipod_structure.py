#!/usr/bin/env python3
"""
Bipod 구조 수정: 강력한 에너지 최소화로 원자 충돌 해결
"""

import sys
from openmm import unit, LangevinIntegrator, Platform, Context, LocalEnergyMinimizer
from openmm import app
from openmm.app import CharmmPsfFile, CharmmCrdFile, CharmmParameterSet

# 파일 경로
BASE_DIR = "/home/pjho3tr/projects/Drug/2026-01-20_Glycogate_Bipod/data/solution builder/openmm"
PSF_FILE = f"{BASE_DIR}/step3_input.psf"
CRD_FILE = f"{BASE_DIR}/step3_input.crd"
TOPPAR_STR = f"{BASE_DIR}/toppar.str"

print("=" * 80)
print("Bipod 구조 수정: 강력한 에너지 최소화")
print("=" * 80)
print()

# 1. 파일 로드
print("Step 1: 파일 로드")
print("-" * 80)
psf = CharmmPsfFile(PSF_FILE)
crd = CharmmCrdFile(CRD_FILE)

# toppar.str에서 파라미터 파일 목록 읽기
with open(TOPPAR_STR, 'r') as f:
    param_files = [line.strip() for line in f if line.strip() and not line.startswith('#')]

# 상대 경로를 절대 경로로 변환
param_files = [f"{BASE_DIR}/{pf}" for pf in param_files]

print(f"✅ PSF: {PSF_FILE}")
print(f"✅ CRD: {CRD_FILE}")
print(f"✅ Parameters: {len(param_files)} files")
print()

# 2. 파라미터 로드
print("Step 2: Force field 로드")
print("-" * 80)
params = CharmmParameterSet(*param_files)
print("✅ Force field 로드 완료")
print()

# 3. Box 정보 설정
print("Step 3: Box 정보 설정")
print("-" * 80)
# CHARMM-GUI step3_pbcsetup.str에서 box 크기 읽기
pbc_setup_file = f"{BASE_DIR}/../step3_pbcsetup.str"
box_a = box_b = box_c = 120.0  # 기본값

try:
    with open(pbc_setup_file, 'r') as f:
        for line in f:
            if 'SET A' in line and '=' in line:
                box_a = float(line.split('=')[1].strip())
            elif 'SET B' in line and '=' in line:
                box_b = float(line.split('=')[1].strip())
            elif 'SET C' in line and '=' in line:
                box_c = float(line.split('=')[1].strip())
    print(f"✅ CHARMM-GUI box 정보 읽음: {pbc_setup_file}")
except:
    print(f"⚠️  {pbc_setup_file} 읽기 실패, 기본값 120 Å 사용")

# PSF에 box 설정
psf.setBox(box_a*unit.angstrom, box_b*unit.angstrom, box_c*unit.angstrom)
print(f"✅ Box 설정: {box_a} x {box_b} x {box_c} Å")
print()

# 4. System 생성
print("Step 4: System 생성")
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

# 5. 초기 에너지 확인
print("Step 5: 초기 에너지 확인")
print("-" * 80)
integrator = LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
platform = Platform.getPlatformByName('CUDA')
context = Context(system, integrator, platform)
context.setPositions(crd.positions)

initial_energy = context.getState(getEnergy=True).getPotentialEnergy()
print(f"초기 에너지: {initial_energy}")
print()

# 6. 단계적 에너지 최소화
print("Step 6: 단계적 에너지 최소화")
print("-" * 80)

# 5-1. 매우 작은 스텝으로 시작 (충돌이 심할 때)
print("Phase 1: 매우 작은 스텝 (tolerance=1000 kJ/mol/nm)")
LocalEnergyMinimizer.minimize(context, tolerance=1000*unit.kilojoule/(unit.mole*unit.nanometer), maxIterations=5000)
energy1 = context.getState(getEnergy=True).getPotentialEnergy()
print(f"  에너지: {energy1}")

# 5-2. 중간 스텝
print("Phase 2: 중간 스텝 (tolerance=100 kJ/mol/nm)")
LocalEnergyMinimizer.minimize(context, tolerance=100*unit.kilojoule/(unit.mole*unit.nanometer), maxIterations=10000)
energy2 = context.getState(getEnergy=True).getPotentialEnergy()
print(f"  에너지: {energy2}")

# 5-3. 정밀 최소화
print("Phase 3: 정밀 최소화 (tolerance=10 kJ/mol/nm)")
LocalEnergyMinimizer.minimize(context, tolerance=10*unit.kilojoule/(unit.mole*unit.nanometer), maxIterations=20000)
final_energy = context.getState(getEnergy=True).getPotentialEnergy()
print(f"  최종 에너지: {final_energy}")
print()

# 7. 결과 저장
print("Step 7: 결과 저장")
print("-" * 80)
state = context.getState(getPositions=True)
positions = state.getPositions()

# PDB 저장
with open(f"{BASE_DIR}/step3_input_minimized.pdb", 'w') as f:
    app.PDBFile.writeFile(psf.topology, positions, f)

# CRD 저장 (OpenMM은 CRD 쓰기를 직접 지원하지 않으므로 PDB만 저장)
print(f"✅ 저장 완료: {BASE_DIR}/step3_input_minimized.pdb")
print()

# 7. 요약
print("=" * 80)
print("요약")
print("=" * 80)
print(f"초기 에너지: {initial_energy}")
print(f"최종 에너지: {final_energy}")
print(f"에너지 감소: {initial_energy - final_energy}")
print()

if final_energy.value_in_unit(unit.kilojoule/unit.mole) < 1e10:
    print("✅ 에너지 최소화 성공!")
    print("   다음 단계: step3_input_minimized.pdb를 사용하여 MD 실행")
else:
    print("⚠️  에너지가 여전히 높습니다.")
    print("   추가 조치 필요 (방법 2 또는 3 참조)")
