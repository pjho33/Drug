#!/bin/bash
#
# OpenMM MD 시뮬레이션 실행 (간단 버전)
#
# CHARMM-GUI가 제공하는 openmm_run.py 사용
#

set -e

echo "================================================================================"
echo "OpenMM MD 시뮬레이션 실행"
echo "================================================================================"
echo

# 경로 설정
WORK_DIR="/home/pjho3/projects/Drug/2026-01-18_Glycogate"
OPENMM_DIR="${WORK_DIR}/data/solution builder/openmm"
RESULTS_DIR="${WORK_DIR}/results/md_1arm_openmm"

# Replica 번호 (명령행 인자)
REPLICA=${1:-1}

echo "작업 디렉토리: ${WORK_DIR}"
echo "OpenMM 파일: ${OPENMM_DIR}"
echo "결과 저장: ${RESULTS_DIR}"
echo "Replica: ${REPLICA}"
echo

# 결과 디렉토리 생성
mkdir -p "${RESULTS_DIR}"

# OpenMM 디렉토리로 이동
cd "${OPENMM_DIR}"

echo "Step 1: OpenMM 환경 확인"
echo "--------------------------------------------------------------------------------"

# OpenMM 버전 확인
python -c "import openmm; print(f'OpenMM version: {openmm.version.version}')" || {
    echo "❌ OpenMM이 설치되지 않았습니다."
    echo "   설치: conda install -c conda-forge openmm"
    exit 1
}

# Platform 확인
python -c "
import openmm
print('사용 가능한 Platform:')
for i in range(openmm.Platform.getNumPlatforms()):
    print(f'  - {openmm.Platform.getPlatform(i).getName()}')
"

echo
echo "✅ OpenMM 환경 확인 완료"
echo

# Step 2: 에너지 최소화 및 평형화
echo "Step 2: 에너지 최소화 및 평형화 (125 ps)"
echo "--------------------------------------------------------------------------------"

INIT="step3_input"
EQUI_PREFIX="step4_equilibration_rep${REPLICA}"

INPUT_PARAM="-t toppar.str -p ${INIT}.psf -c ${INIT}.crd -b sysinfo.dat"

echo "실행 중..."
python -u openmm_run.py -i step4_equilibration.inp ${INPUT_PARAM} \
    -orst "${RESULTS_DIR}/${EQUI_PREFIX}.rst" \
    -odcd "${RESULTS_DIR}/${EQUI_PREFIX}.dcd" \
    > "${RESULTS_DIR}/${EQUI_PREFIX}.out" 2>&1

if [ $? -eq 0 ]; then
    echo "✅ 평형화 완료"
    echo "   출력: ${RESULTS_DIR}/${EQUI_PREFIX}.*"
else
    echo "❌ 평형화 실패"
    echo "   로그: ${RESULTS_DIR}/${EQUI_PREFIX}.out"
    exit 1
fi

echo

# Step 3: Production MD
echo "Step 3: Production MD"
echo "--------------------------------------------------------------------------------"

PROD_PREFIX="step5_production"
PROD_STEP="step5"

# Production 스텝 수 계산
# 기본: 250,000 steps × 4 fs = 1 ns
# 200 ns를 위해서는 200번 반복 필요
# 하지만 일단 10번 (10 ns)으로 테스트

CNT=1
CNTMAX=10  # 10 ns (테스트용)
# CNTMAX=200  # 200 ns (실제 실행용)

echo "Production MD: ${CNTMAX} × 1 ns = ${CNTMAX} ns"
echo

while [ ${CNT} -le ${CNTMAX} ]; do
    PCNT=$((CNT - 1))
    ISTEP="${PROD_STEP}_${CNT}_rep${REPLICA}"
    
    if [ ${CNT} -eq 1 ]; then
        PSTEP="${EQUI_PREFIX}"
    else
        PSTEP="${PROD_STEP}_${PCNT}_rep${REPLICA}"
    fi
    
    echo "  Step ${CNT}/${CNTMAX}: ${ISTEP}"
    
    INPUT_PARAM="-t toppar.str -p ${INIT}.psf -c ${INIT}.crd -irst ${RESULTS_DIR}/${PSTEP}.rst"
    
    python -u openmm_run.py -i ${PROD_PREFIX}.inp ${INPUT_PARAM} \
        -orst "${RESULTS_DIR}/${ISTEP}.rst" \
        -odcd "${RESULTS_DIR}/${ISTEP}.dcd" \
        > "${RESULTS_DIR}/${ISTEP}.out" 2>&1
    
    if [ $? -ne 0 ]; then
        echo "❌ Step ${CNT} 실패"
        echo "   로그: ${RESULTS_DIR}/${ISTEP}.out"
        exit 1
    fi
    
    CNT=$((CNT + 1))
done

echo
echo "✅ Production MD 완료"
echo

# Step 4: Trajectory 병합 (선택사항)
echo "Step 4: Trajectory 파일 목록"
echo "--------------------------------------------------------------------------------"

ls -lh "${RESULTS_DIR}"/*.dcd

echo
echo "================================================================================"
echo "✅ 모든 작업 완료!"
echo "================================================================================"
echo
echo "출력 파일: ${RESULTS_DIR}/"
echo
echo "다음 단계:"
echo "  1. 로그 확인: cat ${RESULTS_DIR}/${EQUI_PREFIX}.out"
echo "  2. Trajectory 분석: MDTraj 또는 MDAnalysis 사용"
echo "  3. 다른 replica 실행: bash $0 2"
echo
