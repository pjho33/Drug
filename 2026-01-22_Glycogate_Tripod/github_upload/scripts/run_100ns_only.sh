#!/usr/bin/env bash
set -euo pipefail

source ~/miniforge3/etc/profile.d/conda.sh
conda activate drug-md

BASE_DIR="/home/pjho3tr/projects/Drug/2026-01-22_Glycogate_Tripod"
OPENMM_DIR="$BASE_DIR/data/solution_builder/openmm"
RESULTS_DIR="$BASE_DIR/results/md_tripod_100ns"

mkdir -p "$RESULTS_DIR"
cd "$OPENMM_DIR"

echo "=========================================="
echo "Tripod 100ns MD 시뮬레이션 시작"
echo "=========================================="

python openmm_run.py \
  --platform CUDA \
  -i step5_production_100ns.inp \
  -p step3_input.psf \
  -c step3_input.crd \
  -irst "$BASE_DIR/results/equilibration/equilibrated.rst" \
  -t toppar.str \
  -odcd "$RESULTS_DIR/tripod_100ns.dcd" \
  -orst "$RESULTS_DIR/tripod_100ns.rst" \
  -opdb "$RESULTS_DIR/tripod_100ns_final.pdb"

echo "=========================================="
echo "✅ 100ns 완료!"
echo "=========================================="
