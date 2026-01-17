#!/bin/bash
# ante-MMPBSA.py 실행 (radii 포함 topology)

MMPBSA_DIR="/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/mmpbsa_amber"

cd "$MMPBSA_DIR"

echo "================================================================================"
echo "ante-MMPBSA.py 실행 (radii 포함 topology)"
echo "================================================================================"
echo ""

echo "입력 파일:"
echo "  step5_input_radii.parm7 (radii 포함)"
echo ""

echo "출력 파일:"
echo "  complex.prmtop"
echo "  receptor.prmtop"
echo "  ligand.prmtop"
echo ""

echo "Masks:"
echo "  Ligand (strip): :306"
echo "  Receptor: :1-305"
echo ""

echo "================================================================================"
echo "실행 중..."
echo "================================================================================"
echo ""

ante-MMPBSA.py \
  -p step5_input_radii.parm7 \
  -c complex.prmtop \
  -r receptor.prmtop \
  -l ligand.prmtop \
  -s ':306' \
  -n ':1-305'

EXIT_CODE=$?

echo ""
echo "================================================================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "✅ ante-MMPBSA.py 성공!"
    echo "================================================================================"
    echo ""
    
    echo "생성된 파일:"
    ls -lh complex.prmtop receptor.prmtop ligand.prmtop 2>/dev/null
    echo ""
    
    echo "다음 단계: MMPBSA.py 실행"
else
    echo "❌ ante-MMPBSA.py 실패 (exit code: $EXIT_CODE)"
    echo "================================================================================"
fi
