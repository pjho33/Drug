#!/usr/bin/env bash
set -euo pipefail

# ✅ 네 파일명에 맞게 여기만 확인/수정
PSF="step3_input.psf"
PDB="step3_input.pdb"
TOPPAR_DIR="toppar"
OPENMM_RUN="openmm_run.py"

run_one () {
  local INP="$1"
  local OUTDIR="$2"

  mkdir -p "$OUTDIR"
  cp -f "$INP" "$OUTDIR/"
  # 입력파일들도 필요하면 복사(선택)
  # cp -f "$PSF" "$PDB" "$OUTDIR/"  # 보통 같은 폴더에서 참조해도 됨

  echo "============================================================"
  echo "RUN: $INP -> $OUTDIR"
  echo "============================================================"

  pushd "$OUTDIR" >/dev/null

  # 원본 폴더의 파일을 그대로 참조(심볼릭 링크)
  ln -sf "../$PSF" .
  ln -sf "../$PDB" .
  ln -sf "../$OPENMM_RUN" .
  ln -sf "../$TOPPAR_DIR" .

  # 실행
  python "$OPENMM_RUN" \
    -i "$(basename "$INP")" \
    -p "$PSF" \
    -c "$PDB" \
    -t "$TOPPAR_DIR" \
    2>&1 | tee run.log

  popd >/dev/null
}

# 1) 10 ns
run_one "../step5_production_10ns.inp" "run_10ns"

echo
echo "✅ 10 ns finished. Check run_10ns/run.log and trajectory outputs."
echo

# 2) 100 ns (10 ns 확인 후 수동으로 실행하는 걸 권장)
# 주석 해제해서 바로 이어서 돌릴 수도 있음.
# run_one "../step5_production_100ns.inp" "run_100ns"
