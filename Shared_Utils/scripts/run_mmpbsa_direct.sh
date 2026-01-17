#!/bin/bash
# Direct MMPBSA.py execution with Amber files

AMBER_DIR="/home/pjho3/projects/Drug/final_complex/GLUT1SDGComplex260110/amber"
MMPBSA_DIR="/home/pjho3/projects/Drug/final_complex/mmpbsa_amber"

mkdir -p "$MMPBSA_DIR"
cd "$MMPBSA_DIR"

echo "================================================================================"
echo "MMPBSA.py Calculation with Amber Format Files"
echo "================================================================================"
echo ""
echo "Input files:"
echo "  Topology: $AMBER_DIR/step5_input.parm7"
echo "  Coordinates: $AMBER_DIR/step5_input.rst7"
echo ""

# Create MMPBSA input file
cat > mmpbsa.in << 'EOF'
Input file for MMPBSA calculation
&general
startframe=1
endframe=1
interval=1
verbose=2
/

&gb
igb=5, saltcon=0.150
/

&pb
istrng=0.150, fillratio=4.0
/
EOF

echo "✅ MMPBSA input file created"
echo ""

# Determine receptor and ligand masks
# Protein is typically residues 16-410 (adjust based on your system)
# Ligand SDG is residue 1
RECEPTOR_MASK=":16-410"
LIGAND_MASK=":1"

echo "Atom masks:"
echo "  Receptor: $RECEPTOR_MASK"
echo "  Ligand: $LIGAND_MASK"
echo ""

echo "================================================================================"
echo "Running MMPBSA.py..."
echo "================================================================================"
echo ""

# Run MMPBSA.py
MMPBSA.py -O \
  -i mmpbsa.in \
  -o FINAL_RESULTS_MMPBSA.dat \
  -sp "$AMBER_DIR/step5_input.parm7" \
  -cp "$AMBER_DIR/step5_input.parm7" \
  -rp "$AMBER_DIR/step5_input.parm7" \
  -lp "$AMBER_DIR/step5_input.parm7" \
  -y "$AMBER_DIR/step5_input.rst7" \
  -srp "$RECEPTOR_MASK" \
  -slp "$LIGAND_MASK" \
  > mmpbsa_stdout.log 2> mmpbsa_stderr.log

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
    echo "Error log (last 50 lines):"
    tail -50 mmpbsa_stderr.log
fi

echo ""
echo "================================================================================"
echo "Output files saved to: $MMPBSA_DIR"
echo "  - FINAL_RESULTS_MMPBSA.dat"
echo "  - mmpbsa_stdout.log"
echo "  - mmpbsa_stderr.log"
echo "================================================================================"
