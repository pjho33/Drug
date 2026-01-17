#!/bin/bash
# 백그라운드 실행 래퍼 스크립트
# 사용법: ./run_bg.sh <script_name.py> [args...]
#
# 기능:
# 1. Drug-MD 환경에서 자동 실행
# 2. 백그라운드 실행
# 3. 실시간 로그 모니터링 가능

if [ $# -eq 0 ]; then
    echo "사용법: $0 <script_name.py> [args...]"
    echo ""
    echo "예시:"
    echo "  $0 verify_radii.py"
    echo "  $0 analyze_trajectory.py --input traj.dcd"
    exit 1
fi

SCRIPT_PATH="$1"
shift  # 첫 번째 인자 제거, 나머지는 스크립트 인자

# 스크립트 이름 추출
SCRIPT_NAME=$(basename "$SCRIPT_PATH" .py)
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# 로그 디렉토리
LOG_DIR="/home/pjho3/projects/Drug/logs"
mkdir -p "$LOG_DIR"

# 로그 파일
LOG_FILE="$LOG_DIR/${SCRIPT_NAME}_${TIMESTAMP}.log"
PID_FILE="$LOG_DIR/${SCRIPT_NAME}_${TIMESTAMP}.pid"

echo "================================================================================"
echo "백그라운드 실행 시작"
echo "================================================================================"
echo "스크립트: $SCRIPT_PATH"
echo "환경: Drug-MD"
echo "로그: $LOG_FILE"
echo "================================================================================"
echo ""

# 백그라운드 실행
nohup conda run -n Drug-MD python "$SCRIPT_PATH" "$@" > "$LOG_FILE" 2>&1 &
BG_PID=$!

# PID 저장
echo $BG_PID > "$PID_FILE"

echo "✅ 백그라운드 실행 시작됨"
echo "   PID: $BG_PID"
echo "   PID 파일: $PID_FILE"
echo ""

# 잠시 대기 후 초기 출력 표시
sleep 2

if ps -p $BG_PID > /dev/null; then
    echo "📊 초기 로그 (최근 20줄):"
    echo "--------------------------------------------------------------------------------"
    tail -20 "$LOG_FILE"
    echo "--------------------------------------------------------------------------------"
    echo ""
    echo "실시간 모니터링:"
    echo "  tail -f $LOG_FILE"
    echo ""
    echo "프로세스 확인:"
    echo "  ps -p $BG_PID"
    echo ""
    echo "프로세스 종료:"
    echo "  kill $BG_PID"
else
    echo "❌ 프로세스가 즉시 종료되었습니다. 로그를 확인하세요:"
    echo "  cat $LOG_FILE"
    exit 1
fi
