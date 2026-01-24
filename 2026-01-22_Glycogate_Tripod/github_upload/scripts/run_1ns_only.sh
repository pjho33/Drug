#!/usr/bin/env bash
set -euo pipefail

source ~/miniforge3/etc/profile.d/conda.sh
conda activate drug-md

BASE_DIR="/home/pjho3tr/projects/Drug/2026-01-22_Glycogate_Tripod"
OPENMM_DIR="$BASE_DIR/data/solution_builder/openmm"
RESULTS_DIR="$BASE_DIR/results/md_tripod_1ns"

mkdir -p "$RESULTS_DIR"
cd "$OPENMM_DIR"

echo "=========================================="
echo "Tripod 1ns MD 시뮬레이션 시작"
echo "=========================================="

python openmm_run.py \
  --platform CUDA \
  -i step5_production_1ns.inp \
  -p step3_input.psf \
  -c step3_input.crd \
  -irst "$BASE_DIR/results/equilibration/equilibrated.rst" \
  -t toppar.str \
  -odcd "$RESULTS_DIR/tripod_1ns.dcd" \
  -orst "$RESULTS_DIR/tripod_1ns.rst" \
  -opdb "$RESULTS_DIR/tripod_1ns_final.pdb"

echo "=========================================="
echo "✅ 1ns 완료!"
echo "=========================================="
