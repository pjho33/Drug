# Phase 3 시뮬레이션 실행 준비 완료

## ✅ 완료된 작업

1. **toppar.str 수정** - Tripod parameters 추가
2. **PDB 복사** - Tripod 포함 PDB 배치
3. **MD 스크립트 작성** - GPU 0, 1 분리

## 🚀 시뮬레이션 실행

```bash
cd /home/pjho3tr/projects/Drug/phase3_with_tripod
bash run_both_simulations.sh
```

## 📊 모니터링

```bash
# 로그 확인
tail -f experimental_gpu0.log
tail -f control_gpu1.log

# 프로세스 확인
ps aux | grep python | grep phase3
```

## 📁 출력 파일

**실험군** (experimental/):
- prod_experimental.dcd - 궤적
- prod_experimental.log - 로그
- prod_experimental.chk - 체크포인트
- prod_experimental_final.pdb - 최종 구조

**대조군** (control/):
- prod_control.dcd - 궤적
- prod_control.log - 로그
- prod_control.chk - 체크포인트
- prod_control_final.pdb - 최종 구조

## ⏱️ 예상 시간

100 ns @ 2 fs timestep = 50,000,000 steps
예상: 24-48시간 (GPU 성능에 따라)
