#!/bin/bash
# 최근 실행 로그 확인 도구

LOG_DIR="/home/pjho3/projects/Drug/logs"

if [ ! -d "$LOG_DIR" ]; then
    echo "로그 디렉토리가 없습니다: $LOG_DIR"
    exit 1
fi

echo "================================================================================"
echo "최근 실행 로그"
echo "================================================================================"
echo ""

# 최근 로그 파일 목록
echo "📁 최근 로그 파일 (최신 10개):"
ls -lt "$LOG_DIR"/*.log 2>/dev/null | head -10 | awk '{print "  " $9 " (" $6" "$7" "$8")"}'
echo ""

# 실행 중인 프로세스
echo "🔄 실행 중인 프로세스:"
RUNNING=0
for pid_file in "$LOG_DIR"/*.pid; do
    if [ -f "$pid_file" ]; then
        PID=$(cat "$pid_file")
        if ps -p $PID > /dev/null 2>&1; then
            SCRIPT_NAME=$(basename "$pid_file" .pid)
            echo "  ✅ $SCRIPT_NAME (PID: $PID)"
            RUNNING=$((RUNNING + 1))
        fi
    fi
done

if [ $RUNNING -eq 0 ]; then
    echo "  (없음)"
fi
echo ""

# 최근 로그 내용
LATEST_LOG=$(ls -t "$LOG_DIR"/*.log 2>/dev/null | head -1)
if [ -n "$LATEST_LOG" ]; then
    echo "📊 최근 로그 내용 ($(basename "$LATEST_LOG")):"
    echo "--------------------------------------------------------------------------------"
    tail -30 "$LATEST_LOG"
    echo "--------------------------------------------------------------------------------"
    echo ""
    echo "전체 로그 보기:"
    echo "  cat $LATEST_LOG"
    echo ""
    echo "실시간 모니터링:"
    echo "  tail -f $LATEST_LOG"
fi
