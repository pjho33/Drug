#!/usr/bin/env python3
"""
OpenMM MD 시뮬레이션 진행 상황 모니터링

실시간으로 로그 파일을 읽어서 진행 상황을 표시합니다.
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
        'total_steps': 0,
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
        
        # 마지막 몇 줄에서 정보 추출
        for line in reversed(lines[-50:]):
            # Progress 정보
            # 예: #"Progress: 10.0%, Time: 100.0 ps, Speed: 150.0 ns/day"
            if 'Progress' in line or 'Step' in line:
                # Progress 추출
                prog_match = re.search(r'Progress[:\s]+([\d.]+)%', line)
                if prog_match:
                    info['progress'] = float(prog_match.group(1))
                
                # Time 추출
                time_match = re.search(r'Time[:\s]+([\d.]+)\s*ps', line)
                if time_match:
                    info['time_ps'] = float(time_match.group(1))
                
                # Speed 추출
                speed_match = re.search(r'Speed[:\s]+([\d.]+)\s*ns/day', line)
                if speed_match:
                    info['speed'] = float(speed_match.group(1))
                
                # Step 추출
                step_match = re.search(r'Step[:\s]+(\d+)', line)
                if step_match:
                    info['step'] = int(step_match.group(1))
            
            # 에너지 정보
            if 'Potential Energy' in line or 'Total Energy' in line:
                energy_match = re.search(r'([-\d.]+)\s*kJ/mol', line)
                if energy_match:
                    info['energy'] = float(energy_match.group(1))
            
            # 온도 정보
            if 'Temperature' in line:
                temp_match = re.search(r'([\d.]+)\s*K', line)
                if temp_match:
                    info['temperature'] = float(temp_match.group(1))
        
        # 파일 수정 시간
        info['last_update'] = datetime.fromtimestamp(os.path.getmtime(log_file))
        
    except Exception as e:
        print(f"로그 파일 파싱 오류: {e}")
    
    return info


def format_time(seconds):
    """초를 시:분:초 형식으로 변환"""
    return str(timedelta(seconds=int(seconds)))


def monitor_simulation(results_dir, replica=1, interval=10):
    """시뮬레이션 모니터링"""
    
    results_path = Path(results_dir)
    
    print("=" * 80)
    print("OpenMM MD 시뮬레이션 모니터링")
    print("=" * 80)
    print()
    print(f"결과 디렉토리: {results_dir}")
    print(f"Replica: {replica}")
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
            print(f"OpenMM MD 모니터링 - Replica {replica}")
            print(f"업데이트: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
            print("=" * 80)
            print()
            
            # 평형화 상태
            equi_log = results_path / f"step4_equilibration_rep{replica}.out"
            print("📊 평형화 (Equilibration)")
            print("-" * 80)
            
            if equi_log.exists():
                equi_info = parse_log_file(str(equi_log))
                if equi_info and equi_info['progress'] > 0:
                    print(f"  ✅ 진행: {equi_info['progress']:.1f}%")
                    print(f"  ⏱️  시간: {equi_info['time_ps']:.1f} ps")
                    if equi_info['speed'] > 0:
                        print(f"  🚀 속도: {equi_info['speed']:.1f} ns/day")
                    if equi_info['temperature']:
                        print(f"  🌡️  온도: {equi_info['temperature']:.1f} K")
                    if equi_info['energy']:
                        print(f"  ⚡ 에너지: {equi_info['energy']:.1f} kJ/mol")
                    
                    if equi_info['progress'] >= 99.9:
                        print("  ✅ 완료!")
                else:
                    print("  ⏳ 시작 중...")
            else:
                print("  ⏳ 대기 중...")
            
            print()
            
            # Production 상태
            print("📊 Production MD")
            print("-" * 80)
            
            # Production 파일 찾기
            prod_files = sorted(results_path.glob(f"step5_*_rep{replica}.out"))
            
            if prod_files:
                total_steps = len(prod_files)
                completed_steps = 0
                current_step = 0
                
                for prod_file in prod_files:
                    step_num = int(re.search(r'step5_(\d+)_rep', prod_file.name).group(1))
                    prod_info = parse_log_file(str(prod_file))
                    
                    if prod_info and prod_info['progress'] >= 99.9:
                        completed_steps += 1
                    elif prod_info and prod_info['progress'] > 0:
                        current_step = step_num
                
                print(f"  📈 완료: {completed_steps} / {total_steps} steps")
                print(f"  ⏱️  총 시간: {completed_steps} ns (목표: {total_steps} ns)")
                
                if current_step > 0:
                    current_log = results_path / f"step5_{current_step}_rep{replica}.out"
                    current_info = parse_log_file(str(current_log))
                    
                    if current_info:
                        print(f"  🔄 현재 Step: {current_step}")
                        print(f"  📊 진행: {current_info['progress']:.1f}%")
                        if current_info['speed'] > 0:
                            print(f"  🚀 속도: {current_info['speed']:.1f} ns/day")
                        if current_info['temperature']:
                            print(f"  🌡️  온도: {current_info['temperature']:.1f} K")
                        
                        # 예상 완료 시간
                        if current_info['speed'] > 0:
                            remaining_ns = total_steps - completed_steps
                            remaining_hours = (remaining_ns / current_info['speed']) * 24
                            eta = datetime.now() + timedelta(hours=remaining_hours)
                            print(f"  ⏰ 예상 완료: {eta.strftime('%Y-%m-%d %H:%M')}")
                            print(f"  ⏳ 남은 시간: {format_time(remaining_hours * 3600)}")
                
                # 진행률 바
                progress_pct = (completed_steps / total_steps) * 100 if total_steps > 0 else 0
                bar_length = 50
                filled = int(bar_length * progress_pct / 100)
                bar = '█' * filled + '░' * (bar_length - filled)
                print(f"  [{bar}] {progress_pct:.1f}%")
                
            else:
                print("  ⏳ 시작 대기 중...")
            
            print()
            
            # 디스크 사용량
            print("💾 디스크 사용량")
            print("-" * 80)
            
            dcd_files = list(results_path.glob(f"*_rep{replica}.dcd"))
            total_size = sum(f.stat().st_size for f in dcd_files) / (1024**3)  # GB
            
            print(f"  📁 DCD 파일: {len(dcd_files)}개")
            print(f"  💿 총 크기: {total_size:.2f} GB")
            
            print()
            print("=" * 80)
            print(f"실행 시간: {format_time(time.time() - start_time)}")
            print("Ctrl+C를 눌러 종료")
            print("=" * 80)
            
            time.sleep(interval)
            
    except KeyboardInterrupt:
        print("\n\n모니터링 종료")


def main():
    """메인 함수"""
    
    # 기본 경로
    results_dir = "/home/pjho3/projects/Drug/2026-01-18_Glycogate/results/md_1arm_openmm"
    
    # 명령행 인자
    replica = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    interval = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    
    monitor_simulation(results_dir, replica, interval)


if __name__ == "__main__":
    main()
