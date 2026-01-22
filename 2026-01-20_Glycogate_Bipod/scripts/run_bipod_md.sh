#!/usr/bin/env bash
set -euo pipefail

# Bipod MD 시뮬레이션 실행 스크립트
# 최소화된 구조를 사용하여 10ns MD 실행

source ~/miniforge3/etc/profile.d/conda.sh
conda activate drug-md

OPENMM_DIR="/home/pjho3tr/projects/Drug/2026-01-20_Glycogate_Bipod/data/solution builder/openmm"
RESULTS_DIR="/home/pjho3tr/projects/Drug/2026-01-20_Glycogate_Bipod/results/md_bipod_10ns"

mkdir -p "$RESULTS_DIR"

cd "$OPENMM_DIR"

echo "=========================================="
echo "Bipod MD 시뮬레이션 시작"
echo "=========================================="
echo "입력 구조: step3_input_minimized.pdb"
echo "시뮬레이션 길이: 10 ns"
echo "출력 디렉토리: $RESULTS_DIR"
echo ""

# 최소화된 PDB를 CRD로 변환 (OpenMM은 PDB를 직접 사용 가능)
# 하지만 openmm_run.py는 CRD를 요구하므로, 원본 CRD를 그대로 사용
# (최소화는 이미 완료되었고, topology는 동일)

nohup python openmm_run.py \
  -i step5_production_10ns.inp \
  -p step3_input.psf \
  -c step3_input.crd \
  -t toppar.str \
  -odcd "$RESULTS_DIR/bipod_10ns.dcd" \
  -orst "$RESULTS_DIR/bipod_10ns.rst" \
  -opdb "$RESULTS_DIR/bipod_10ns_final.pdb" \
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
