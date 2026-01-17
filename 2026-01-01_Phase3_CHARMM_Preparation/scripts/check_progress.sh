#!/bin/bash
echo "=========================================="
echo "Phase 3 MD Simulation Progress"
echo "=========================================="
echo ""

# Check if simulations are running
echo "=== Running Processes ==="
ps aux | grep "run_experimental_gpu0.py\|run_control_gpu1.py" | grep -v grep

echo ""
echo "=== GPU Status ==="
nvidia-smi --query-gpu=index,name,utilization.gpu,memory.used,memory.total --format=csv,noheader

echo ""
echo "=== Experimental (GPU 0) Progress ==="
if [ -f experimental/prod_experimental.log ]; then
    LAST_LINE=$(tail -1 experimental/prod_experimental.log)
    echo "Latest: $LAST_LINE"
    
    # Calculate progress
    if [ -f experimental/prod_experimental.log ]; then
        CURRENT_STEP=$(tail -100 experimental/prod_experimental.log | grep -v "^#" | grep "," | tail -1 | cut -d',' -f1)
        if [ ! -z "$CURRENT_STEP" ]; then
            TOTAL_STEPS=50000000
            PROGRESS=$(echo "scale=2; $CURRENT_STEP / $TOTAL_STEPS * 100" | bc)
            TIME_NS=$(echo "scale=2; $CURRENT_STEP * 0.002 / 1000" | bc)
            echo "Progress: ${PROGRESS}% (${TIME_NS} ns / 100 ns)"
        fi
    fi
else
    echo "No log file yet"
fi

echo ""
echo "=== Control (GPU 1) Progress ==="
if [ -f control/prod_control.log ]; then
    LAST_LINE=$(tail -1 control/prod_control.log)
    echo "Latest: $LAST_LINE"
    
    # Calculate progress
    if [ -f control/prod_control.log ]; then
        CURRENT_STEP=$(tail -100 control/prod_control.log | grep -v "^#" | grep "," | tail -1 | cut -d',' -f1)
        if [ ! -z "$CURRENT_STEP" ]; then
            TOTAL_STEPS=50000000
            PROGRESS=$(echo "scale=2; $CURRENT_STEP / $TOTAL_STEPS * 100" | bc)
            TIME_NS=$(echo "scale=2; $CURRENT_STEP * 0.002 / 1000" | bc)
            echo "Progress: ${PROGRESS}% (${TIME_NS} ns / 100 ns)"
        fi
    fi
else
    echo "No log file yet"
fi

echo ""
echo "=== Checkpoint Files ==="
ls -lh experimental/*.chk control/*.chk 2>/dev/null || echo "No checkpoints yet"

echo ""
echo "=== DCD Files ==="
ls -lh experimental/*.dcd control/*.dcd 2>/dev/null || echo "No trajectories yet"

echo ""
echo "=========================================="
