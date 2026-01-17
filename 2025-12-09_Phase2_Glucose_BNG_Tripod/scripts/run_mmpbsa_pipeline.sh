#!/bin/bash
#
# Complete MM-GBSA Pipeline for CHARMM-GUI/OpenMM outputs
# Safe route: PSF → dry PSF → Amber → MM-GBSA
#

set -e  # Exit on error

echo "======================================================================"
echo "MM-GBSA Pipeline for Glycosylated GLUT1-SDG Complex"
echo "======================================================================"

# Configuration
TOP_PSF="/home/pjho3tr/Downloads/GlycatedGLUT1SDGComplex260107/openmm/step5_input.psf"
TOP_PDB="/home/pjho3tr/Downloads/GlycatedGLUT1SDGComplex260107/openmm/step5_input.pdb"
TOPPAR_DIR="/home/pjho3tr/Downloads/GlycatedGLUT1SDGComplex260107/toppar"
LIGAND_TOPPAR="/home/pjho3tr/Downloads/GlycatedGLUT1SDGComplex260107/sdg"
TRAJ_DCD="/home/pjho3tr/Downloads/GlycatedGLUT1SDGComplex260107/openmm/production.dcd"

LIG_RESNAME="SDG"
GLYCAN_RESNAMES="NAG BMA MAN GAL FUC SIA NEU5AC NDG BGLC A2G A2M"

WORK_DIR="/home/pjho3tr/projects/Drug/mmpbsa"
SCRIPT_DIR="/home/pjho3tr/projects/Drug/scripts"

cd "$WORK_DIR"

echo ""
echo "======================================================================"
echo "Step 1: Create dry trajectory (protein + glycan + SDG)"
echo "======================================================================"
if [ ! -f "dry_complex.dcd" ]; then
    echo "⏭️  Using existing dry trajectory from previous step..."
    if [ ! -f "dry_complex.dcd" ]; then
        echo "❌ Error: dry_complex.dcd not found. Run make_dry_traj.py first."
        exit 1
    fi
else
    echo "✅ Dry trajectory already exists"
fi

echo ""
echo "======================================================================"
echo "Step 2: Create dry PSF (exclude water/ions/lipids)"
echo "======================================================================"
if [ ! -f "dry_complex.psf" ]; then
    conda run -n drug-md python "$SCRIPT_DIR/make_dry_psf.py" \
        --psf "$TOP_PSF" \
        --pdb "$TOP_PDB" \
        --out "dry_complex" \
        --lig "$LIG_RESNAME" \
        --glycan $GLYCAN_RESNAMES
    echo "✅ Dry PSF created"
else
    echo "⏭️  Dry PSF already exists, skipping..."
fi

echo ""
echo "======================================================================"
echo "Step 3: Convert dry PSF to Amber prmtop"
echo "======================================================================"
if [ ! -f "complex.parm7" ]; then
    conda run -n drug-md python "$SCRIPT_DIR/psf_to_amber.py" \
        --psf "dry_complex.psf" \
        --pdb "dry_complex.pdb" \
        --toppar "$TOPPAR_DIR" \
        --ligand-toppar "$LIGAND_TOPPAR" \
        --out "complex"
    echo "✅ Amber topology created"
else
    echo "⏭️  Amber topology already exists, skipping..."
fi

echo ""
echo "======================================================================"
echo "Step 4: Add mbondi2 radii for GB calculations"
echo "======================================================================"
if [ ! -f "complex_radii.parm7" ]; then
    conda run -n drug-md python "$SCRIPT_DIR/add_radii.py" \
        --in "complex.parm7" \
        --rst "complex.rst7" \
        --out "complex_radii"
    echo "✅ Radii added"
else
    echo "⏭️  Radii file already exists, skipping..."
fi

echo ""
echo "======================================================================"
echo "Step 5: Run ante-MMPBSA (split receptor/ligand)"
echo "======================================================================"
if [ ! -f "receptor.prmtop" ]; then
    conda run -n drug-md ante-MMPBSA.py \
        -p complex_radii.parm7 \
        -c complex_radii.prmtop \
        -r receptor.prmtop \
        -l ligand.prmtop \
        -s ":$LIG_RESNAME" \
        -n ":$LIG_RESNAME"
    echo "✅ Receptor/ligand topologies created"
else
    echo "⏭️  Receptor/ligand topologies already exist, skipping..."
fi

echo ""
echo "======================================================================"
echo "Step 6: Run MM-GBSA calculation"
echo "======================================================================"
echo "Running MM-GBSA on full trajectory (10,000 frames)..."
echo "This will take 2-4 hours. Running in background..."

nohup conda run -n drug-md MMPBSA.py -O \
    -i "$SCRIPT_DIR/mmpbsa_gb.in" \
    -cp complex_radii.prmtop \
    -rp receptor.prmtop \
    -lp ligand.prmtop \
    -y dry_complex.dcd \
    -o FINAL_GBSA.dat \
    -eo perframe_GBSA.csv > mmpbsa_run.log 2>&1 &

MMPBSA_PID=$!
echo "✅ MM-GBSA started (PID: $MMPBSA_PID)"
echo "   Monitor progress: tail -f $WORK_DIR/mmpbsa_run.log"

echo ""
echo "======================================================================"
echo "Pipeline Setup Complete!"
echo "======================================================================"
echo ""
echo "Next steps:"
echo "  1. Wait for MM-GBSA to complete (check mmpbsa_run.log)"
echo "  2. Analyze results in FINAL_GBSA.dat"
echo "  3. Optional: Split trajectory into 0-25ns and 40-100ns for comparison"
echo ""
