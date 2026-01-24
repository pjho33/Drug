#!/usr/bin/env bash
set -euo pipefail

source ~/miniforge3/etc/profile.d/conda.sh
conda activate drug-md

OPENMM_DIR="/home/pjho3tr/projects/Drug/2026-01-22_Glycogate_Tripod/data/solution_builder/openmm"
RESULTS_DIR="/home/pjho3tr/projects/Drug/2026-01-22_Glycogate_Tripod/results/md_tripod_100ns"

mkdir -p "$RESULTS_DIR"
cd "$OPENMM_DIR"

echo "=========================================="
echo "Tripod MD 시뮬레이션 시작 (100ns)"
echo "=========================================="

nohup python openmm_run.py \
  -i step5_production_100ns.inp \
  -p step3_input.psf \
  -c step3_input.crd \
  -t toppar.str \
  -odcd "$RESULTS_DIR/tripod_100ns.dcd" \
  -orst "$RESULTS_DIR/tripod_100ns.rst" \
  -opdb "$RESULTS_DIR/tripod_100ns_final.pdb" \
  > "$RESULTS_DIR/run.log" 2>&1 &

PID=$!
echo "✅ 100ns MD 시작 (PID: $PID)"
echo "$PID" > "$RESULTS_DIR/md.pid"
