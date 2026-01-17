#!/bin/bash
# Monitor MD simulation progress

echo "=========================================="
echo "MD Simulation Monitor"
echo "=========================================="
echo ""

# Check if simulation is running
if pgrep -f "run_control_md.py" > /dev/null; then
    echo "✅ Simulation is RUNNING"
    echo ""
    
    # Show process info
    echo "Process info:"
    ps aux | grep "run_control_md.py" | grep -v grep
    echo ""
else
    echo "⚠️  Simulation is NOT running"
    echo ""
fi

# Show last 30 lines of log
echo "=========================================="
echo "Latest log output (last 30 lines):"
echo "=========================================="
if [ -f simulation.log ]; then
    tail -n 30 simulation.log
else
    echo "Log file not found yet..."
fi

echo ""
echo "=========================================="
echo "Output files:"
echo "=========================================="
ls -lh /home/pjho3/projects/Drug/control_md_simulation/results/ 2>/dev/null || echo "Results directory not created yet..."

echo ""
echo "=========================================="
echo "Commands:"
echo "  Watch log: tail -f simulation.log"
echo "  Check process: ps aux | grep run_control_md"
echo "  Stop simulation: pkill -f run_control_md.py"
echo "=========================================="
