#!/bin/bash
# Monitor control complex MD simulation status

WORK_DIR="/home/pjho3/projects/Drug/final_complex/controlcomplex/openmm"
LOG_FILE="/home/pjho3/projects/Drug/final_complex/controlcomplex/simulation.log"

echo "=========================================="
echo "Control Complex MD Simulation Status"
echo "=========================================="
echo ""

# Check if simulation is running
if pgrep -f "run_control_simulation.sh" > /dev/null; then
    echo "✅ Simulation is RUNNING"
    echo ""
else
    echo "⚠️  Simulation process not found"
    echo ""
fi

# Check log file
if [ -f "$LOG_FILE" ]; then
    echo "📋 Latest log entries:"
    echo "------------------------------------------"
    tail -20 "$LOG_FILE"
    echo ""
else
    echo "⚠️  Log file not found: $LOG_FILE"
    echo ""
fi

# Check equilibration progress
echo "=========================================="
echo "Equilibration Progress"
echo "=========================================="
cd "$WORK_DIR"

for i in {1..6}; do
    step="step6.${i}_equilibration"
    if [ -f "${step}.rst" ]; then
        echo "✅ ${step} - COMPLETED"
        if [ -f "${step}.out" ]; then
            # Check for errors
            if grep -q "ERROR" "${step}.out"; then
                echo "   ⚠️  Errors found in output"
            fi
        fi
    elif [ -f "${step}.out" ]; then
        echo "🔄 ${step} - IN PROGRESS"
        # Show last few lines
        echo "   Last output:"
        tail -3 "${step}.out" | sed 's/^/   /'
    else
        echo "⏳ ${step} - PENDING"
    fi
done

echo ""

# Check production progress
echo "=========================================="
echo "Production Progress"
echo "=========================================="

if [ -f "step7_1.rst" ]; then
    echo "✅ step7_1 (Production) - COMPLETED"
elif [ -f "step7_1.out" ]; then
    echo "🔄 step7_1 (Production) - IN PROGRESS"
    echo ""
    echo "Last output:"
    tail -5 "step7_1.out" | sed 's/^/   /'
else
    echo "⏳ step7_1 (Production) - PENDING"
fi

echo ""

# Check output files
echo "=========================================="
echo "Output Files"
echo "=========================================="
echo ""
echo "DCD trajectory files:"
ls -lh step*.dcd 2>/dev/null | awk '{print "  " $9 " (" $5 ")"}'

echo ""
echo "Restart files:"
ls -lh step*.rst 2>/dev/null | awk '{print "  " $9 " (" $5 ")"}'

echo ""
echo "=========================================="
echo "To view full log: tail -f $LOG_FILE"
echo "To check again: bash $0"
echo "=========================================="
