#!/bin/bash
#
# 1-Arm PEG24-Glc MD 시스템 준비
#
# CHARMM-GUI 파일을 사용하여 GROMACS MD 시뮬레이션 준비
#

set -e

echo "================================================================================"
echo "1-Arm MD 시스템 준비"
echo "================================================================================"
echo

# 경로 설정
WORK_DIR="/home/pjho3/projects/Drug/2026-01-18_Glycogate"
DATA_DIR="${WORK_DIR}/data"
CHARMM_GUI_DIR="${DATA_DIR}/charmm-gui-6871763698"
GROMACS_DIR="${CHARMM_GUI_DIR}/gromacs"
RESULTS_DIR="${WORK_DIR}/results/md_1arm"

# 결과 디렉토리 생성
mkdir -p "${RESULTS_DIR}"

echo "작업 디렉토리: ${WORK_DIR}"
echo "CHARMM-GUI 출력: ${CHARMM_GUI_DIR}"
echo "결과 저장: ${RESULTS_DIR}"
echo

# CHARMM-GUI 파일 확인
echo "Step 1: CHARMM-GUI 파일 확인"
echo "--------------------------------------------------------------------------------"

if [ ! -f "${GROMACS_DIR}/topol.top" ]; then
    echo "❌ topol.top 파일이 없습니다."
    echo "   CHARMM-GUI Solution Builder로 시스템을 구축해주세요."
    exit 1
fi

if [ ! -f "${GROMACS_DIR}/LIG.itp" ]; then
    echo "❌ LIG.itp 파일이 없습니다."
    exit 1
fi

echo "✅ Topology 파일 확인 완료"
echo

# 구조 파일 확인
echo "Step 2: 구조 파일 확인"
echo "--------------------------------------------------------------------------------"

if [ -f "${CHARMM_GUI_DIR}/ligandrm.pdb" ]; then
    echo "✅ ligandrm.pdb 발견"
    STRUCT_FILE="${CHARMM_GUI_DIR}/ligandrm.pdb"
elif [ -f "${DATA_DIR}/1arm_peg24_glc.pdb" ]; then
    echo "✅ 1arm_peg24_glc.pdb 사용"
    STRUCT_FILE="${DATA_DIR}/1arm_peg24_glc.pdb"
else
    echo "❌ 구조 파일이 없습니다."
    exit 1
fi

echo "구조 파일: ${STRUCT_FILE}"
echo

# 다음 단계 안내
echo "================================================================================"
echo "✅ 파일 확인 완료!"
echo "================================================================================"
echo
echo "다음 단계:"
echo
echo "1. CHARMM-GUI Solution Builder로 수용액 시스템 구축 (권장)"
echo "   URL: http://www.charmm-gui.org/?doc=input/solution"
echo
echo "   설정:"
echo "   - PDB 파일: ${STRUCT_FILE}"
echo "   - Water: TIP3P"
echo "   - Ion: 150 mM NaCl"
echo "   - Box: 12 nm cubic"
echo "   - Output: GROMACS"
echo
echo "2. 또는 수동으로 시스템 구축:"
echo "   cd ${RESULTS_DIR}"
echo "   gmx editconf -f ${STRUCT_FILE} -o box.gro -c -d 6.0 -bt cubic"
echo "   gmx solvate -cp box.gro -cs spc216.gro -o solvated.gro -p topol.top"
echo "   # ... (이온 추가, 에너지 최소화 등)"
echo
echo "3. MD 시뮬레이션 실행"
echo "   - 200-500 ns"
echo "   - 3 replica (다른 seed)"
echo

echo "참고 문서: docs/04_MD_Simulation_Setup.md"
echo
