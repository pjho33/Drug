#!/usr/bin/env bash
# Tripod 시뮬레이션 실시간 모니터링

BASE_DIR="/home/pjho3tr/projects/Drug/2026-01-22_Glycogate_Tripod/results"

echo "=========================================="
echo "Tripod MD 시뮬레이션 모니터링"
echo "=========================================="
echo ""

# 프로세스 확인
if ps aux | grep -q "[o]penmm_run.py"; then
    echo "🟢 OpenMM 실행 중"
    ps aux | grep "[o]penmm_run.py" | awk '{print "   PID:", $2, "CPU:", $3"%", "MEM:", $4"%"}'
else
    echo "⚪ OpenMM 프로세스 없음"
fi
echo ""

# GPU 사용 확인
echo "🎮 GPU 상태:"
nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv,noheader 2>/dev/null || echo "   GPU 사용 중인 프로세스 없음"
echo ""

# 각 단계별 상태
for stage in equilibration md_tripod_1ns md_tripod_10ns md_tripod_100ns; do
    LOG="$BASE_DIR/$stage/run.log"
    if [ -f "$LOG" ]; then
        echo "📊 $stage:"
        
        # 에러 확인
        if grep -q "Error\|Exception\|NaN" "$LOG"; then
            echo "   ❌ 에러 발생!"
            grep -i "error\|exception\|nan" "$LOG" | tail -3
        else
            # 마지막 몇 줄 확인
            echo "   최근 로그:"
            tail -3 "$LOG" | sed 's/^/   /'
        fi
        echo ""
    fi
done

# Master 로그
echo "📝 Master 로그:"
if [ -f "$BASE_DIR/master.log" ]; then
    tail -5 "$BASE_DIR/master.log" | sed 's/^/   /'
else
    echo "   로그 없음"
fi

echo ""
echo "=========================================="
echo "💡 실시간 로그 확인:"
echo "   tail -f $BASE_DIR/master.log"
echo "   tail -f $BASE_DIR/equilibration/run.log"
echo "   tail -f $BASE_DIR/md_tripod_1ns/run.log"
echo "=========================================="
