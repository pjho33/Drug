#!/bin/bash

# Arguments: <Target_PDB> <Ligand_SMILES> <Project_Name>
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_receptor_pdb> <ligand_smiles_string> <project_name>"
    exit 1
fi

RECEPTOR_PDB=$1
LIGAND_SMILES=$2
PROJECT_NAME=$3

OUTPUT_DIR="results/$PROJECT_NAME"
CLEANED_RECEPTOR="$OUTPUT_DIR/cleaned_receptor.pdb"
DOCKING_OUTPUT_DIR="$OUTPUT_DIR/docking_output/complex_0"

# --- Setup ---
echo "📂 Creating Output Directories: $OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$DOCKING_OUTPUT_DIR"

echo "=========================================================="
echo "🚀 [Start] Auto Drug Discovery Pipeline (MANUAL): $PROJECT_NAME"
echo "   Target: $RECEPTOR_PDB"
echo "   Drug: $LIGAND_SMILES"
echo "=========================================================="

# ----------------------------------------------------------
# 🔹 [Step 1] Preparing Receptor (Environment: dd_final)
# ----------------------------------------------------------
echo ""
echo "🔹 [Step 1] Preparing Receptor (Environment: dd_final)..."
conda run -n dd_final python scripts/01_prepare_receptor.py "$RECEPTOR_PDB" "$CLEANED_RECEPTOR"
if [ $? -ne 0 ]; then echo "❌ Step 1 Failed!"; exit 1; fi

# ----------------------------------------------------------
# 🔹 [Step 2] Positioning Tripod Manually
# ----------------------------------------------------------
echo ""
echo "🔹 [Step 2] Positioning Tripod Manually (AI Skipped)..."
conda run -n dd_final python scripts/02_manual_docking.py "$CLEANED_RECEPTOR" "$LIGAND_SMILES" "$DOCKING_OUTPUT_DIR"
if [ ! -f "$DOCKING_OUTPUT_DIR/rank1.sdf" ]; then echo "❌ Step 2 Failed! Manual placement file not found."; exit 1; fi


# --- Determine Ligand Input for MD/Scoring ---
LIGAND_INPUT="$DOCKING_OUTPUT_DIR/rank1.sdf"
TRAJECTORY_FILE="$OUTPUT_DIR/validation_trajectory.pdb"

# ----------------------------------------------------------
# 🔹 [Step 3] Running MD Validation (Environment: dd_final)
# **스크립트 파일명 03_run_simulation2.py로 변경**
# ----------------------------------------------------------
echo ""
echo "🔹 [Step 3] Running MD Validation (Environment: dd_final)..."
conda run -n dd_final python scripts/03_run_simulation2.py "$CLEANED_RECEPTOR" "$LIGAND_INPUT" "$TRAJECTORY_FILE"
if [ $? -ne 0 ]; then echo "❌ Step 3 Failed!"; exit 1; fi

# ----------------------------------------------------------
# 🔹 [Step 4] Scoring Binding Energy (MM-GBSA)
# ----------------------------------------------------------
echo ""
echo "🔹 [Step 4] Scoring Binding Energy (MM-GBSA)..."
conda run -n dd_final python scripts/04_calc_energy.py "$TRAJECTORY_FILE" "$LIGAND_INPUT" "$OUTPUT_DIR/binding_score.csv"
if [ $? -ne 0 ]; then echo "❌ Step 4 Failed!"; exit 1; fi

# --- Success ---
echo "=========================================================="
echo "🎉 [Success] Pipeline Complete! (Manual Mode)"
echo "   📂 Trajectory: $TRAJECTORY_FILE"
echo "   📊 Final Score: $OUTPUT_DIR/binding_score.csv"
echo "=========================================================="