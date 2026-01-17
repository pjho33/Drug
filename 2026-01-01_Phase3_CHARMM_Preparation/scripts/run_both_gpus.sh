#!/bin/bash

echo "=========================================="
echo "Starting Phase 3: Dual GPU Simulation"
echo "=========================================="
echo "GPU 0: Glycosylated system"
echo "GPU 1: Control system"
echo "=========================================="

# Activate conda environment
source /home/pjho3tr/miniforge3/etc/profile.d/conda.sh
conda activate drug-md

cd /home/pjho3tr/projects/Drug/phase3_glycosylation

# Start both simulations
python run_glyco_gpu0.py > glyco_gpu0.log 2>&1 &
PID1=$!

python run_control_gpu1.py > control_gpu1.log 2>&1 &
PID2=$!

echo "Started Glycosylated on GPU 0 (PID: $PID1)"
echo "Started Control on GPU 1 (PID: $PID2)"
echo ""
echo "✓ 두 시뮬레이션 시작됨"
echo ""
echo "진행 상황 확인:"
echo "  bash check_progress.sh"
echo ""
echo "로그 확인:"
echo "  tail -f glyco_gpu0.log"
echo "  tail -f control_gpu1.log"
echo ""
echo "중지:"
echo "  kill $PID1 $PID2"
