#!/bin/bash
#
# Create dry prmtop matching dry_complex.dcd (7432 atoms)
# Using cpptraj to strip from original CHARMM topology
#

set -e

echo "======================================================================"
echo "Creating dry prmtop from CHARMM PSF (matching 7432 atoms DCD)"
echo "======================================================================"

# Paths
ORIG_DIR="/home/pjho3tr/Downloads/GlycatedGLUT1SDGComplex260107/openmm"
ORIG_PSF="$ORIG_DIR/step5_input.psf"
ORIG_PDB="$ORIG_DIR/step5_input.pdb"
TOPPAR_DIR="/home/pjho3tr/Downloads/GlycatedGLUT1SDGComplex260107/toppar"
SDG_DIR="/home/pjho3tr/Downloads/GlycatedGLUT1SDGComplex260107/sdg"

WORK_DIR="/home/pjho3tr/projects/Drug/mmpbsa"
cd "$WORK_DIR"

echo ""
echo "Step 1: Convert CHARMM PSF to Amber prmtop (full system)"
echo "======================================================================"

# Use the existing conversion script
conda run -n drug-md python /home/pjho3tr/projects/Drug/scripts/psf_to_amber.py \
    --psf "$ORIG_PSF" \
    --pdb "$ORIG_PDB" \
    --toppar "$TOPPAR_DIR" \
    --ligand-toppar "$SDG_DIR" \
    --out full_system

echo ""
echo "Step 2: Strip to dry system (protein + glycan + SDG)"
echo "======================================================================"

# Strip water, ions, lipids using cpptraj
# Keep: protein + SDG + glycans (NAG, BMA, MAN, GAL, FUC, SIA, etc.)
cat > strip_to_dry.in <<'EOF'
# Strip everything except protein, glycans, and SDG
strip :WAT,TIP3,HOH
strip :NA,CL,K,CA,MG,POT,CLA,SOD,CAL
strip :POPC,POPE,POPS,POPG,CHL1,CHOL
strip :DLPC,DPPC,DOPC,DOPE,DOPS,DPPE,DPPG

# Write dry topology
parmwrite out dry_complex.parm7
EOF

conda run -n drug-md cpptraj -p full_system.parm7 -i strip_to_dry.in

echo ""
echo "Step 3: Verify atom count"
echo "======================================================================"

conda run -n drug-md python -c "
import parmed as pmd
parm = pmd.load_file('dry_complex.parm7')
print(f'Dry prmtop atoms: {len(parm.atoms):,}')
print(f'Expected (from DCD): 7,432')
if len(parm.atoms) == 7432:
    print('✅ Perfect match!')
else:
    print(f'⚠️  Mismatch: {len(parm.atoms)} vs 7432')
"

echo ""
echo "======================================================================"
echo "Complete!"
echo "======================================================================"
