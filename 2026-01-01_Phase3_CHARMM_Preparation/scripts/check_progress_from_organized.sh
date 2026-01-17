#!/bin/bash
# Phase 3 진행 상황 모니터링

echo "=========================================="
echo "Phase 3 시뮬레이션 진행 상황"
echo "=========================================="
echo ""

# GPU 프로세스 확인
echo "=== GPU 프로세스 ==="
ps aux | grep -E "run_glyco|run_control" | grep -v grep || echo "실행 중인 프로세스 없음"
echo ""

# 실험모델 상태
echo "=== 실험모델 (Glycosylated) - GPU 0 ==="
if [ -f /home/pjho3tr/projects/Drug/phase3_glycosylation/glycosylated_new_final/prod_glyco_tripod_log.csv ]; then
    LINES=$(wc -l < /home/pjho3tr/projects/Drug/phase3_glycosylation/glycosylated_new_final/prod_glyco_tripod_log.csv)
    if [ $LINES -gt 1 ]; then
        LAST=$(tail -1 /home/pjho3tr/projects/Drug/phase3_glycosylation/glycosylated_new_final/prod_glyco_tripod_log.csv)
        STEP=$(echo $LAST | cut -d',' -f1)
        TIME=$(echo $LAST | cut -d',' -f2)
        TEMP=$(echo $LAST | cut -d',' -f4)
        SPEED=$(echo $LAST | cut -d',' -f5)
        
        PROGRESS=$(echo "scale=2; $STEP / 50000000 * 100" | bc)
        NS=$(echo "scale=1; $TIME / 1000" | bc)
        
        echo "  진행률: ${PROGRESS}% (${NS} ns / 100 ns)"
        echo "  Step: $STEP / 50,000,000"
        echo "  온도: ${TEMP} K"
        echo "  속도: ${SPEED} ns/day"
        
        # 남은 시간 계산
        if [ ! -z "$SPEED" ] && [ "$SPEED" != "0" ]; then
            REMAINING_NS=$(echo "100 - $NS" | bc)
            REMAINING_DAYS=$(echo "scale=2; $REMAINING_NS / $SPEED" | bc)
            echo "  예상 완료: ${REMAINING_DAYS} 일"
        fi
    else
        echo "  초기화 중..."
    fi
else
    echo "  아직 시작 안 됨"
fi
echo ""

# 대조모델 상태
echo "=== 대조모델 (Control) - GPU 1 ==="
if [ -f /home/pjho3tr/projects/Drug/phase3_glycosylation/control_final/prod_control_tripod_log.csv ]; then
    LINES=$(wc -l < /home/pjho3tr/projects/Drug/phase3_glycosylation/control_final/prod_control_tripod_log.csv)
    if [ $LINES -gt 1 ]; then
        LAST=$(tail -1 /home/pjho3tr/projects/Drug/phase3_glycosylation/control_final/prod_control_tripod_log.csv)
        STEP=$(echo $LAST | cut -d',' -f1)
        TIME=$(echo $LAST | cut -d',' -f2)
        TEMP=$(echo $LAST | cut -d',' -f4)
        SPEED=$(echo $LAST | cut -d',' -f5)
        
        PROGRESS=$(echo "scale=2; $STEP / 50000000 * 100" | bc)
        NS=$(echo "scale=1; $TIME / 1000" | bc)
        
        echo "  진행률: ${PROGRESS}% (${NS} ns / 100 ns)"
        echo "  Step: $STEP / 50,000,000"
        echo "  온도: ${TEMP} K"
        echo "  속도: ${SPEED} ns/day"
        
        if [ ! -z "$SPEED" ] && [ "$SPEED" != "0" ]; then
            REMAINING_NS=$(echo "100 - $NS" | bc)
            REMAINING_DAYS=$(echo "scale=2; $REMAINING_NS / $SPEED" | bc)
            echo "  예상 완료: ${REMAINING_DAYS} 일"
        fi
    else
        echo "  초기화 중..."
    fi
else
    echo "  아직 시작 안 됨"
fi
echo ""

# 파일 크기
echo "=== 파일 크기 ==="
ls -lh /home/pjho3tr/projects/Drug/phase3_glycosylation/*/prod_*.dcd 2>/dev/null | awk '{print $9, $5}' || echo "DCD 파일 없음"
echo ""

echo "=========================================="
echo "사용법: bash check_progress.sh"
echo "=========================================="
