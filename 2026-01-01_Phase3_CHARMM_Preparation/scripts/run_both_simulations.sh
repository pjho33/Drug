#!/bin/bash
echo "Phase 3: Dual GPU MD Simulation"
echo "Experimental (GPU 0) | Control (GPU 1)"

source /home/pjho3tr/miniforge3/etc/profile.d/conda.sh
conda activate drug-md

cd /home/pjho3tr/projects/Drug/phase3_with_tripod

python run_experimental_gpu0.py > experimental_gpu0.log 2>&1 &
PID1=$!

python run_control_gpu1.py > control_gpu1.log 2>&1 &
PID2=$!

echo "Started: PID $PID1 (GPU 0) | PID $PID2 (GPU 1)"
echo "Monitor: tail -f experimental_gpu0.log"
echo "Stop: kill $PID1 $PID2"
