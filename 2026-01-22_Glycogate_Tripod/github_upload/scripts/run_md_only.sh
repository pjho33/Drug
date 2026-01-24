#!/usr/bin/env bash
set -euo pipefail

# Tripod: Equilibration 완료 후 MD만 실행 (1ns → 10ns → 100ns)

source ~/miniforge3/etc/profile.d/conda.sh
conda activate drug-md

BASE_DIR="/home/pjho3tr/projects/Drug/2026-01-22_Glycogate_Tripod"
OPENMM_DIR="$BASE_DIR/data/solution_builder/openmm"
RESULTS_BASE="$BASE_DIR/results"

cd "$OPENMM_DIR"

# 1ns 시뮬레이션
echo "=========================================="
echo "Tripod 1ns MD 시뮬레이션 시작"
echo "=========================================="

RESULTS_DIR="$RESULTS_BASE/md_tripod_1ns"
mkdir -p "$RESULTS_DIR"

python openmm_run.py \
  --platform CUDA \
  -i step5_production_1ns.inp \
  -p step3_input.psf \
  -c step3_input.crd \
  -irst "$RESULTS_BASE/equilibration/equilibrated.rst" \
  -t toppar.str \
  -odcd "$RESULTS_DIR/tripod_1ns.dcd" \
  -orst "$RESULTS_DIR/tripod_1ns.rst" \
  -opdb "$RESULTS_DIR/tripod_1ns_final.pdb" \
  > "$RESULTS_DIR/run.log" 2>&1

echo "✅ 1ns 완료"
echo ""

# 10ns 시뮬레이션
echo "=========================================="
echo "Tripod 10ns MD 시뮬레이션 시작"
echo "=========================================="

RESULTS_DIR="$RESULTS_BASE/md_tripod_10ns"
mkdir -p "$RESULTS_DIR"

python openmm_run.py \
  --platform CUDA \
  -i step5_production_10ns.inp \
  -p step3_input.psf \
  -c step3_input.crd \
  -irst "$RESULTS_BASE/equilibration/equilibrated.rst" \
  -t toppar.str \
  -odcd "$RESULTS_DIR/tripod_10ns.dcd" \
  -orst "$RESULTS_DIR/tripod_10ns.rst" \
  -opdb "$RESULTS_DIR/tripod_10ns_final.pdb" \
  > "$RESULTS_DIR/run.log" 2>&1

echo "✅ 10ns 완료"
echo ""

# 100ns 시뮬레이션
echo "=========================================="
echo "Tripod 100ns MD 시뮬레이션 시작"
echo "=========================================="

RESULTS_DIR="$RESULTS_BASE/md_tripod_100ns"
mkdir -p "$RESULTS_DIR"

python openmm_run.py \
  --platform CUDA \
  -i step5_production_100ns.inp \
  -p step3_input.psf \
  -c step3_input.crd \
  -irst "$RESULTS_BASE/equilibration/equilibrated.rst" \
  -t toppar.str \
  -odcd "$RESULTS_DIR/tripod_100ns.dcd" \
  -orst "$RESULTS_DIR/tripod_100ns.rst" \
  -opdb "$RESULTS_DIR/tripod_100ns_final.pdb" \
  > "$RESULTS_DIR/run.log" 2>&1

echo "✅ 100ns 완료"
echo "=========================================="
echo "모든 Tripod 시뮬레이션 완료!"
echo "=========================================="
