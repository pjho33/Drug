#!/bin/bash

# ========================================================
# [사용법] ./run_pipeline.sh <PDB파일> <SMILES> <프로젝트명>
# Version 4.0: MM-GBSA Ligand Parameter Fix
# ========================================================

# ... (앞부분 동일) ...

# 1. 입력 파라미터 확인
INPUT_PDB=$1
INPUT_SMILES=$2
PROJECT_NAME=$3

if [ -z "$PROJECT_NAME" ]; then
    echo "❌ Usage: ./run_pipeline.sh <pdb_file> <smiles> <project_name>"
    exit 1
fi

# 2. 폴더 생성
BASE_DIR=$(pwd)
RESULTS_DIR="$BASE_DIR/results/$PROJECT_NAME"
LOG_DIR="$RESULTS_DIR/logs"

echo "📂 Creating Output Directories: $RESULTS_DIR"
mkdir -p "$RESULTS_DIR"
mkdir -p "$LOG_DIR"

LOG_FILE="$LOG_DIR/pipeline.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=========================================================="
echo "🚀 [Start] Auto Drug Discovery Pipeline: $PROJECT_NAME"
echo "   Target: $INPUT_PDB"
echo "   Drug: $INPUT_SMILES"
echo "=========================================================="

# --------------------------------------------------------
# [Step 1] 수용체 준비 (dd_core)
# --------------------------------------------------------
echo ""
echo "🔹 [Step 1] Preparing Receptor (Environment: dd_core)..."
eval "$(conda shell.bash hook)"
conda activate dd_core

CLEANED_PDB="$RESULTS_DIR/cleaned_receptor.pdb"
python scripts/01_prepare_receptor.py "$INPUT_PDB" "$CLEANED_PDB"

if [ ! -f "$CLEANED_PDB" ]; then
    echo "❌ Step 1 Failed!"
    exit 1
fi

# --------------------------------------------------------
# [Step 2] AI 도킹 (dd_ai)
# --------------------------------------------------------
echo ""
echo "🔹 [Step 2] Running DiffDock (Environment: dd_ai)..."
conda activate dd_ai

DOCKING_OUT_DIR="$RESULTS_DIR/docking_output"
python scripts/02_run_docking.py "$CLEANED_PDB" "$INPUT_SMILES" "$DOCKING_OUT_DIR"

FOUND_FILE=$(find "$DOCKING_OUT_DIR" -name "rank1*.sdf" | head -n 1)

if [ -z "$FOUND_FILE" ]; then
    echo "❌ Step 2 Failed! Docking result not found."
    exit 1
else
    echo "   ✅ Found candidate: $FOUND_FILE"
    BEST_LIGAND="$RESULTS_DIR/best_docking_pose.sdf"
    cp "$FOUND_FILE" "$BEST_LIGAND"
fi

# --------------------------------------------------------
# [Step 3] MD 검증 (dd_core)
# --------------------------------------------------------
echo ""
echo "🔹 [Step 3] Running Validation (Environment: dd_core)..."
conda activate dd_core

VALIDATION_TRAJ="$RESULTS_DIR/validation_trajectory.pdb"
python scripts/03_run_simulation.py "$CLEANED_PDB" "$BEST_LIGAND" "$VALIDATION_TRAJ"

if [ ! -f "$VALIDATION_TRAJ" ]; then
    echo "❌ Step 3 Failed!"
    exit 1
fi

# --------------------------------------------------------
# [Step 4] 결합 에너지 계산 (Fix: Pass Ligand File)
# --------------------------------------------------------
echo ""
echo "🔹 [Step 4] Scoring Binding Energy (MM-GBSA)..."

SCORE_FILE="$RESULTS_DIR/binding_score.csv"
# ✅ 수정된 부분: BEST_LIGAND를 인자로 넘겨줍니다!
python scripts/04_calc_energy.py "$VALIDATION_TRAJ" "$BEST_LIGAND" "$SCORE_FILE"

echo ""
echo "=========================================================="
echo "🎉 [Success] Pipeline Complete!"
echo "   📂 Trajectory: $VALIDATION_TRAJ"
echo "   📊 Final Score: $SCORE_FILE"
echo "=========================================================="
tail -n 2 "$SCORE_FILE"