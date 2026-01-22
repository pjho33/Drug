#!/usr/bin/env bash
set -euo pipefail

# Bipod 100ns MD 시뮬레이션 실행 스크립트

source ~/miniforge3/etc/profile.d/conda.sh
conda activate drug-md

OPENMM_DIR="/home/pjho3tr/projects/Drug/2026-01-20_Glycogate_Bipod/data/solution builder/openmm"
RESULTS_DIR="/home/pjho3tr/projects/Drug/2026-01-20_Glycogate_Bipod/results/md_bipod_100ns"

mkdir -p "$RESULTS_DIR"

cd "$OPENMM_DIR"

echo "=========================================="
echo "Bipod MD 시뮬레이션 시작 (100ns)"
echo "=========================================="
echo "입력 구조: step3_input.crd"
echo "시뮬레이션 길이: 100 ns"
echo "출력 디렉토리: $RESULTS_DIR"
echo ""

nohup python openmm_run.py \
  -i step5_production_100ns.inp \
  -p step3_input.psf \
  -c step3_input.crd \
  -t toppar.str \
  -odcd "$RESULTS_DIR/bipod_100ns.dcd" \
  -orst "$RESULTS_DIR/bipod_100ns.rst" \
  -opdb "$RESULTS_DIR/bipod_100ns_final.pdb" \
  > "$RESULTS_DIR/run.log" 2>&1 &

PID=$!
echo "✅ MD 시뮬레이션 시작됨 (PID: $PID)"
echo "$PID" > "$RESULTS_DIR/md.pid"
echo ""
echo "진행 상황 확인:"
echo "  tail -f $RESULTS_DIR/run.log"
echo ""
echo "프로세스 종료:"
echo "  kill $PID"
