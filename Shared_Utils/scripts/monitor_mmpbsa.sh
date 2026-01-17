#!/bin/bash
# Monitor MMPBSA.py progress

MMPBSA_DIR="/home/pjho3/projects/Drug/final_complex/mmpbsa_amber"

echo "================================================================================"
echo "MMPBSA.py Progress Monitor"
echo "================================================================================"
echo ""

# Check if MMPBSA is running
MMPBSA_PID=$(pgrep -f "MMPBSA.py")
if [ -n "$MMPBSA_PID" ]; then
    echo "✅ MMPBSA.py is running (PID: $MMPBSA_PID)"
else
    echo "⚠️  MMPBSA.py is not currently running"
fi

echo ""
echo "================================================================================"
echo "Recent output (last 30 lines):"
echo "================================================================================"

if [ -f "$MMPBSA_DIR/mmpbsa_stdout.log" ]; then
    tail -30 "$MMPBSA_DIR/mmpbsa_stdout.log"
else
    echo "No output log found yet"
fi

echo ""
echo "================================================================================"
echo "Errors (if any):"
echo "================================================================================"

if [ -f "$MMPBSA_DIR/mmpbsa_stderr.log" ]; then
    if [ -s "$MMPBSA_DIR/mmpbsa_stderr.log" ]; then
        tail -20 "$MMPBSA_DIR/mmpbsa_stderr.log"
    else
        echo "No errors"
    fi
else
    echo "No error log found yet"
fi

echo ""
echo "================================================================================"
echo "Results status:"
echo "================================================================================"

if [ -f "$MMPBSA_DIR/FINAL_RESULTS_MMPBSA.dat" ]; then
    echo "✅ Results file exists!"
    echo ""
    cat "$MMPBSA_DIR/FINAL_RESULTS_MMPBSA.dat"
else
    echo "⏳ Results not ready yet"
fi

echo ""
echo "================================================================================"
echo "To monitor in real-time, use:"
echo "  tail -f $MMPBSA_DIR/mmpbsa_stdout.log"
echo "================================================================================"
