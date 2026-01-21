#!/bin/bash
#
# CHARMM-GUI OpenMM으로 200ns MD 실행
#

OPENMM_DIR="/home/pjho3tr/projects/Drug/2026-01-18_Glycogate/data/solution builder/openmm"
RESULTS_DIR="/home/pjho3tr/projects/Drug/2026-01-18_Glycogate/results/md_200ns_rep1"

mkdir -p "$RESULTS_DIR"

echo "=========================================="
echo "TRIS-PEG24-Lglucose 200ns MD (Replica 1)"
echo "=========================================="
echo ""

# step5_production.inp를 200ns로 수정
cd "$OPENMM_DIR"
cp step5_production.inp step5_production_200ns.inp

# nstep을 100,000,000으로 변경 (200ns)
sed -i 's/nstep.*=.*/nstep       = 100000000/' step5_production_200ns.inp

echo "✅ 설정 파일 생성: step5_production_200ns.inp"
echo ""

# 배경 실행
# -i: input parameter file
# -p: topology file (PSF) ← 중요!
# -c: coordinate file (CRD)
# -t: toppar stream file ← 중요!
nohup conda run -n drug-md python openmm_run.py \
  -i step5_production_200ns.inp \
  -p step3_input.psf \
  -c step3_input.crd \
  -t toppar.str \
  -odcd "$RESULTS_DIR/md_rep1.dcd" \
  -opdb "$RESULTS_DIR/md_rep1_final.pdb" \
  > "$RESULTS_DIR/run.log" 2>&1 &

PID=$!
echo $PID > "$RESULTS_DIR/simulation.pid"

echo "✅ 시뮬레이션 시작됨"
echo "   PID: $PID"
echo "   작업 디렉토리: $OPENMM_DIR"
echo "   결과 디렉토리: $RESULTS_DIR"
echo "   로그: $RESULTS_DIR/run.log"
echo ""
echo "진행 상황 확인:"
echo "  tail -f $RESULTS_DIR/run.log"
echo ""
echo "중지:"
echo "  kill $PID"
echo ""
