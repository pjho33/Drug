#!/bin/bash
#
# CHARMM-GUI OpenMM으로 200ns MD 실행
#

OPENMM_DIR="/home/pjho3tr/projects/Drug/2026-01-18_Glycogate/data/solution builder/openmm"
RESULTS_DIR="/home/pjho3tr/projects/Drug/2026-01-18_Glycogate/results/md_200ns_rep1"

mkdir -p "$RESULTS_DIR"
cd "$OPENMM_DIR"

echo "=========================================="
echo "TRIS-PEG24-Lglucose 200ns MD (Replica 1)"
echo "=========================================="
echo ""

# step5_production.inp를 200ns로 수정
cat > step5_production_200ns.inp << 'EOF'
# Production MD (200 ns)
nstep = 100000000
dt = 0.002

nstdcd = 5000
nstout = 5000
nstlog = 5000

nbupdate = 10
nstlist = 20

cutoff = 12.0
switchdist = 10.0
pairlistdist = 13.5

pmegrid = 1.0
pmeftol = 1e-5

temp = 300.0
fric_coeff = 1.0

pres = 1.0
barostatfreq = 25

constraints = HBonds
rigidwater = yes
EOF

# 배경 실행
nohup conda run -n drug-md python openmm_run.py \
  -i step5_production_200ns.inp \
  -p toppar.str \
  -c step3_input.crd \
  -t step3_input.psf \
  -o "$RESULTS_DIR/md_rep1" \
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
