#!/bin/bash
#
# MD 시뮬레이션 진행 상황 빠른 확인
#

RESULTS_DIR="/home/pjho3/projects/Drug/2026-01-18_Glycogate/results/md_1arm_openmm"
REPLICA=${1:-1}

echo "================================================================================"
echo "MD 시뮬레이션 진행 상황 (Replica ${REPLICA})"
echo "================================================================================"
echo

# 백그라운드 프로세스 확인
echo "📊 실행 중인 프로세스:"
echo "--------------------------------------------------------------------------------"
ps aux | grep -E "openmm_run.py|09_run_openmm" | grep -v grep || echo "  ⚠️  실행 중인 프로세스 없음"
echo

# 평형화 상태
echo "📊 평형화 (Equilibration):"
echo "--------------------------------------------------------------------------------"
EQUI_LOG="${RESULTS_DIR}/step4_equilibration_rep${REPLICA}.out"
if [ -f "${EQUI_LOG}" ]; then
    echo "  ✅ 로그 파일 존재"
    tail -5 "${EQUI_LOG}" | grep -E "Progress|Step|Time" || echo "  ⏳ 진행 중..."
else
    echo "  ⏳ 아직 시작 안 됨"
fi
echo

# Production 상태
echo "📊 Production MD:"
echo "--------------------------------------------------------------------------------"
PROD_COUNT=$(ls -1 ${RESULTS_DIR}/step5_*_rep${REPLICA}.dcd 2>/dev/null | wc -l)
echo "  📁 완료된 스텝: ${PROD_COUNT}"

if [ ${PROD_COUNT} -gt 0 ]; then
    LATEST_STEP=$(ls -1 ${RESULTS_DIR}/step5_*_rep${REPLICA}.out 2>/dev/null | tail -1)
    if [ -n "${LATEST_STEP}" ]; then
        echo "  📄 최신 로그: $(basename ${LATEST_STEP})"
        tail -3 "${LATEST_STEP}" | grep -E "Progress|Step|Time" || echo "  ⏳ 진행 중..."
    fi
fi
echo

# 디스크 사용량
echo "💾 디스크 사용량:"
echo "--------------------------------------------------------------------------------"
du -sh ${RESULTS_DIR} 2>/dev/null || echo "  디렉토리 없음"
echo

# 최근 로그
echo "📝 최근 로그 (마지막 10줄):"
echo "--------------------------------------------------------------------------------"
MAIN_LOG="${RESULTS_DIR}/simulation_rep${REPLICA}.log"
if [ -f "${MAIN_LOG}" ]; then
    tail -10 "${MAIN_LOG}"
else
    echo "  로그 파일 없음"
fi
echo

echo "================================================================================"
echo "상세 모니터링: python scripts/10_monitor_md.py ${REPLICA}"
echo "================================================================================"
