# 백그라운드 실행 시스템 사용 가이드

모든 Python 스크립트를 **Drug-MD 환경**에서 **백그라운드**로 실행하고 **실시간 모니터링**할 수 있는 시스템입니다.

---

## 🚀 빠른 시작

### 1. 스크립트 백그라운드 실행

```bash
cd /home/pjho3/projects/Drug/scripts
./run_bg.sh <스크립트명.py>
```

**예시:**
```bash
./run_bg.sh verify_radii.py
./run_bg.sh analyze_trajectory.py
./run_bg.sh calculate_mmpbsa.py
```

### 2. 실행 상태 확인

```bash
./check_logs.sh
```

이 명령으로 다음을 확인할 수 있습니다:
- 최근 로그 파일 목록
- 현재 실행 중인 프로세스
- 최근 로그 내용 (30줄)

---

## 📁 로그 관리

### 로그 위치
```
/home/pjho3/projects/Drug/logs/
├── verify_radii_20260111_210000.log
├── verify_radii_20260111_210000.pid
├── analyze_trajectory_20260111_210500.log
└── analyze_trajectory_20260111_210500.pid
```

### 로그 파일 명명 규칙
```
<스크립트명>_<날짜>_<시간>.log
<스크립트명>_<날짜>_<시간>.pid
```

---

## 🔍 모니터링

### 실시간 로그 보기

```bash
# 최근 로그 실시간 모니터링
tail -f /home/pjho3/projects/Drug/logs/<최근로그파일>.log
```

### 특정 스크립트 로그 찾기

```bash
ls -lt /home/pjho3/projects/Drug/logs/verify_radii*.log
```

### 프로세스 확인

```bash
# PID 파일에서 프로세스 ID 확인
cat /home/pjho3/projects/Drug/logs/<스크립트>_<타임스탬프>.pid

# 프로세스 실행 여부 확인
ps -p <PID>
```

---

## 🛑 프로세스 제어

### 프로세스 종료

```bash
# PID로 종료
kill <PID>

# 강제 종료
kill -9 <PID>
```

### 모든 실행 중인 프로세스 확인

```bash
ps aux | grep "conda run -n Drug-MD"
```

---

## 💡 사용 예시

### 시나리오 1: RMSD 분석 실행

```bash
# 1. 백그라운드 실행
./run_bg.sh analyze_rmsd.py

# 출력:
# ✅ 백그라운드 실행 시작됨
#    PID: 12345
#    로그: /home/pjho3/projects/Drug/logs/analyze_rmsd_20260111_210000.log

# 2. 잠시 후 상태 확인
./check_logs.sh

# 3. 실시간 모니터링 (필요시)
tail -f /home/pjho3/projects/Drug/logs/analyze_rmsd_20260111_210000.log
```

### 시나리오 2: 여러 스크립트 동시 실행

```bash
# 여러 분석을 동시에 실행
./run_bg.sh calc_rmsd.py
./run_bg.sh calc_rmsf.py
./run_bg.sh calc_hbonds.py

# 모두 확인
./check_logs.sh
```

### 시나리오 3: 긴 작업 (MMPBSA 등)

```bash
# 1. 실행
./run_bg.sh run_mmpbsa.py

# 2. 주기적으로 확인 (예: 10분마다)
watch -n 600 ./check_logs.sh

# 또는 실시간 모니터링
tail -f /home/pjho3/projects/Drug/logs/run_mmpbsa_*.log
```

---

## ⚙️ 자동 기능

### 1. Drug-MD 환경 자동 활성화
- 모든 스크립트는 자동으로 `Drug-MD` conda 환경에서 실행됩니다
- 환경 활성화를 신경 쓸 필요 없음

### 2. 로그 자동 저장
- 모든 출력(stdout, stderr)이 자동으로 로그 파일에 저장됩니다
- 타임스탬프가 자동으로 추가되어 구분이 쉽습니다

### 3. PID 자동 추적
- 각 실행의 PID가 자동으로 저장됩니다
- 나중에 프로세스를 쉽게 찾고 제어할 수 있습니다

---

## 🔧 고급 사용법

### 인자가 있는 스크립트 실행

```bash
./run_bg.sh analyze_trajectory.py --input traj.dcd --output results.csv
```

### 여러 작업 순차 실행 (스크립트로)

```bash
cat > run_pipeline.sh << 'EOF'
#!/bin/bash
./run_bg.sh step1_prepare.py
sleep 60  # 1분 대기
./run_bg.sh step2_analyze.py
sleep 60
./run_bg.sh step3_visualize.py
EOF

chmod +x run_pipeline.sh
./run_pipeline.sh
```

---

## 📊 로그 정리

### 오래된 로그 삭제

```bash
# 7일 이상 된 로그 삭제
find /home/pjho3/projects/Drug/logs -name "*.log" -mtime +7 -delete
find /home/pjho3/projects/Drug/logs -name "*.pid" -mtime +7 -delete
```

### 특정 스크립트 로그만 삭제

```bash
rm /home/pjho3/projects/Drug/logs/test_*.log
rm /home/pjho3/projects/Drug/logs/test_*.pid
```

---

## 🆘 문제 해결

### Q: 스크립트가 즉시 종료됩니다
**A:** 로그 파일을 확인하세요:
```bash
cat /home/pjho3/projects/Drug/logs/<최근로그>.log
```

### Q: conda 환경이 작동하지 않습니다
**A:** conda 초기화 확인:
```bash
conda info --envs
conda activate Drug-MD  # 수동 테스트
```

### Q: 로그가 업데이트되지 않습니다
**A:** 프로세스가 실행 중인지 확인:
```bash
ps aux | grep python
```

---

## 📝 체크리스트

실행 전:
- [ ] 스크립트가 `/home/pjho3/projects/Drug/scripts/`에 있는가?
- [ ] 스크립트에 실행 권한이 있는가? (`chmod +x`)

실행 중:
- [ ] `./check_logs.sh`로 주기적으로 확인
- [ ] 에러 발생 시 로그 파일 확인

실행 후:
- [ ] 결과 파일 확인
- [ ] 필요 없는 로그 정리

---

## 🎯 핵심 명령어 요약

```bash
# 실행
./run_bg.sh <script.py>

# 상태 확인
./check_logs.sh

# 실시간 모니터링
tail -f /home/pjho3/projects/Drug/logs/<최근로그>.log

# 프로세스 종료
kill <PID>
```

---

**Last Updated:** 2026-01-11
