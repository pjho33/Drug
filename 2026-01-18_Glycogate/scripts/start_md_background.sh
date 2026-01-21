#!/bin/bash
#
# TRIS-PEG24-Lglucose 200ns MD 시뮬레이션 배경 실행
# Replica 1
#

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
RESULTS_DIR="/home/pjho3tr/projects/Drug/2026-01-18_Glycogate/results/md_200ns_rep1"

mkdir -p "$RESULTS_DIR"

echo "=========================================="
echo "TRIS-PEG24-Lglucose 200ns MD (Replica 1)"
echo "=========================================="
echo ""
echo "배경 실행 시작..."
echo ""

# conda 환경에서 nohup으로 배경 실행
nohup conda run -n drug-md python "$SCRIPT_DIR/run_md_200ns_rep1.py" > "$RESULTS_DIR/run.log" 2>&1 &

PID=$!
echo $PID > "$RESULTS_DIR/simulation.pid"

echo "✅ 시뮬레이션 시작됨"
echo "   PID: $PID"
echo "   로그: $RESULTS_DIR/run.log"
echo "   진행 상황: $RESULTS_DIR/md_rep1.log"
echo ""
echo "진행 상황 확인:"
echo "  tail -f $RESULTS_DIR/run.log"
echo "  tail -f $RESULTS_DIR/md_rep1.log"
echo ""
echo "모니터링 스크립트:"
echo "  python3 $SCRIPT_DIR/monitor_md_rep1.py"
echo ""
echo "중지하려면:"
echo "  kill $PID"
echo ""
