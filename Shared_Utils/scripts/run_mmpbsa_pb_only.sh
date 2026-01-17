#!/bin/bash
# Run MMPBSA.py with PB only (no GB)

AMBER_DIR="/home/pjho3/projects/Drug/final_complex/GLUT1SDGComplex260110/amber"
MMPBSA_DIR="/home/pjho3/projects/Drug/final_complex/mmpbsa_amber"

cd "$MMPBSA_DIR"

echo "================================================================================"
echo "MMPBSA.py Calculation (PB only, no GB)"
echo "================================================================================"
echo ""

# Create MMPBSA input file - PB only
cat > mmpbsa_pb.in << 'EOF'
MMPBSA calculation with PB only
&general
startframe=1
endframe=1
interval=1
verbose=2
/

&pb
istrng=0.150, fillratio=4.0
/
EOF

echo "✅ MMPBSA input file created (PB only)"
echo ""

echo "================================================================================"
echo "Running MMPBSA.py..."
echo "================================================================================"
echo ""

# Run MMPBSA.py
MMPBSA.py -O \
  -i mmpbsa_pb.in \
  -o FINAL_RESULTS_MMPBSA.dat \
  -sp "$MMPBSA_DIR/complex.prmtop" \
  -cp "$MMPBSA_DIR/complex.prmtop" \
  -rp "$MMPBSA_DIR/receptor.prmtop" \
  -lp "$MMPBSA_DIR/ligand.prmtop" \
  -y "$AMBER_DIR/step5_input.rst7" \
  > mmpbsa_pb_stdout.log 2> mmpbsa_pb_stderr.log

EXIT_CODE=$?

echo ""
echo "================================================================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "✅ MMPBSA.py completed successfully!"
    echo "================================================================================"
    echo ""
    
    if [ -f "FINAL_RESULTS_MMPBSA.dat" ]; then
        echo "RESULTS:"
        echo "================================================================================"
        cat FINAL_RESULTS_MMPBSA.dat
        echo ""
    fi
else
    echo "❌ MMPBSA.py failed with exit code $EXIT_CODE"
    echo "================================================================================"
    echo ""
    echo "Error log:"
    cat mmpbsa_pb_stderr.log
fi

echo ""
echo "================================================================================"
echo "Output files:"
echo "  - FINAL_RESULTS_MMPBSA.dat"
echo "  - mmpbsa_pb_stdout.log"
echo "  - mmpbsa_pb_stderr.log"
echo "================================================================================"
