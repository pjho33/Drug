#!/usr/bin/env python3
"""
TRIS-PEG24-Lglucose MD 시뮬레이션 진행 상황 모니터링 (Replica 1)
"""

import os
import sys
import time
import re
from pathlib import Path
from datetime import datetime, timedelta

def parse_log_file(log_file):
    """로그 파일에서 진행 상황 파싱"""
    
    if not os.path.exists(log_file):
        return None
    
    info = {
        'step': 0,
        'time_ps': 0.0,
        'speed': 0.0,
        'progress': 0.0,
        'energy': None,
        'temperature': None,
        'last_update': None
    }
    
    try:
        with open(log_file, 'r') as f:
            lines = f.readlines()
        
        if len(lines) < 2:
            return info
        
        # 헤더 찾기
        header_idx = -1
        for i, line in enumerate(lines):
            if '#"Step"' in line or 'Step' in line and 'Time' in line:
                header_idx = i
                break
        
        if header_idx == -1:
            return info
        
        # 마지막 데이터 라인
        for line in reversed(lines[header_idx+1:]):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split(',')
            if len(parts) >= 3:
                try:
                    info['step'] = int(parts[0])
                    info['time_ps'] = float(parts[1])
                    
                    # Speed 찾기 (마지막 컬럼)
                    if len(parts) >= 10:
                        info['speed'] = float(parts[-1])
                    
                    # Temperature
                    if len(parts) >= 7:
                        info['temperature'] = float(parts[6])
                    
                    # Energy
                    if len(parts) >= 5:
                        info['energy'] = float(parts[4])
                    
                    # Progress 계산 (200 ns = 100,000,000 steps)
                    info['progress'] = (info['step'] / 100000000) * 100
                    
                    break
                except (ValueError, IndexError):
                    continue
        
        # 파일 수정 시간
        info['last_update'] = datetime.fromtimestamp(os.path.getmtime(log_file))
        
    except Exception as e:
        print(f"로그 파일 파싱 오류: {e}")
    
    return info


def format_time(seconds):
    """초를 시:분:초 형식으로 변환"""
    return str(timedelta(seconds=int(seconds)))


def monitor_simulation(results_dir, interval=10):
    """시뮬레이션 모니터링"""
    
    results_path = Path(results_dir)
    log_file = results_path / "md_rep1.log"
    run_log = results_path / "run.log"
    pid_file = results_path / "simulation.pid"
    
    print("=" * 80)
    print("TRIS-PEG24-Lglucose 200ns MD 모니터링 (Replica 1)")
    print("=" * 80)
    print()
    print(f"결과 디렉토리: {results_dir}")
    print(f"업데이트 간격: {interval}초")
    print()
    print("Ctrl+C를 눌러 종료")
    print("=" * 80)
    print()
    
    start_time = time.time()
    
    try:
        while True:
            os.system('clear' if os.name == 'posix' else 'cls')
            
            print("=" * 80)
            print(f"TRIS-PEG24-Lglucose 200ns MD - Replica 1")
            print(f"업데이트: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
            print("=" * 80)
            print()
            
            # PID 확인
            if pid_file.exists():
                with open(pid_file, 'r') as f:
                    pid = f.read().strip()
                print(f"🔹 PID: {pid}")
                
                # 프로세스 실행 확인
                try:
                    os.kill(int(pid), 0)
                    print(f"✅ 프로세스 실행 중")
                except OSError:
                    print(f"⚠️  프로세스 종료됨")
            else:
                print("⚠️  PID 파일 없음")
            
            print()
            
            # MD 진행 상황
            print("📊 시뮬레이션 진행 상황")
            print("-" * 80)
            
            if log_file.exists():
                info = parse_log_file(str(log_file))
                
                if info and info['step'] > 0:
                    time_ns = info['time_ps'] / 1000
                    
                    print(f"  ⏱️  시간: {time_ns:.2f} ns / 200.00 ns")
                    print(f"  📈 진행: {info['progress']:.2f}%")
                    print(f"  🔢 Step: {info['step']:,} / 100,000,000")
                    
                    if info['speed'] > 0:
                        print(f"  🚀 속도: {info['speed']:.2f} ns/day")
                        
                        # 예상 완료 시간
                        remaining_ns = 200 - time_ns
                        remaining_days = remaining_ns / info['speed']
                        remaining_hours = remaining_days * 24
                        eta = datetime.now() + timedelta(hours=remaining_hours)
                        
                        print(f"  ⏰ 예상 완료: {eta.strftime('%Y-%m-%d %H:%M')}")
                        print(f"  ⏳ 남은 시간: {format_time(remaining_hours * 3600)}")
                    
                    if info['temperature']:
                        print(f"  🌡️  온도: {info['temperature']:.1f} K")
                    
                    if info['energy']:
                        print(f"  ⚡ 에너지: {info['energy']:.1f} kJ/mol")
                    
                    # 진행률 바
                    bar_length = 50
                    filled = int(bar_length * info['progress'] / 100)
                    bar = '█' * filled + '░' * (bar_length - filled)
                    print(f"  [{bar}] {info['progress']:.1f}%")
                    
                    if info['last_update']:
                        elapsed = (datetime.now() - info['last_update']).total_seconds()
                        if elapsed < 60:
                            print(f"  🕐 마지막 업데이트: {elapsed:.0f}초 전")
                        else:
                            print(f"  🕐 마지막 업데이트: {elapsed/60:.1f}분 전")
                else:
                    print("  ⏳ 시작 중...")
            else:
                print("  ⏳ 로그 파일 대기 중...")
            
            print()
            
            # 디스크 사용량
            print("💾 디스크 사용량")
            print("-" * 80)
            
            dcd_file = results_path / "md_rep1.dcd"
            if dcd_file.exists():
                size_gb = dcd_file.stat().st_size / (1024**3)
                print(f"  📁 DCD: {size_gb:.2f} GB")
            else:
                print("  📁 DCD: 생성 대기 중...")
            
            chk_file = results_path / "md_rep1.chk"
            if chk_file.exists():
                size_mb = chk_file.stat().st_size / (1024**2)
                print(f"  💾 Checkpoint: {size_mb:.2f} MB")
            
            print()
            
            # 로그 파일 정보
            print("📄 로그 파일")
            print("-" * 80)
            print(f"  실행 로그: {run_log}")
            print(f"  MD 로그: {log_file}")
            print()
            print("  확인: tail -f {run_log}")
            print()
            
            print("=" * 80)
            print(f"모니터링 시간: {format_time(time.time() - start_time)}")
            print("Ctrl+C를 눌러 종료")
            print("=" * 80)
            
            time.sleep(interval)
            
    except KeyboardInterrupt:
        print("\n\n모니터링 종료")


def main():
    """메인 함수"""
    
    results_dir = "/home/pjho3tr/projects/Drug/2026-01-18_Glycogate/results/md_200ns_rep1"
    interval = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    
    monitor_simulation(results_dir, interval)


if __name__ == "__main__":
    main()
