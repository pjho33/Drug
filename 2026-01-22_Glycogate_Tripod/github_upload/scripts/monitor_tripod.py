#!/usr/bin/env python3
"""
Tripod MD 시뮬레이션 진행상황 모니터링
"""
import os
import time
import subprocess

BASE_DIR = "/home/pjho3tr/projects/Drug/2026-01-22_Glycogate_Tripod/results"
SIMULATIONS = ["md_tripod_1ns", "md_tripod_10ns", "md_tripod_100ns"]
EXPECTED_FRAMES = {
    "md_tripod_1ns": 10,      # 250k steps / 25k stride = 10 frames
    "md_tripod_10ns": 100,    # 2.5M steps / 25k stride = 100 frames
    "md_tripod_100ns": 1000   # 25M steps / 25k stride = 1000 frames
}

def get_dcd_size(sim_name):
    """DCD 파일 크기 확인"""
    dcd_file = os.path.join(BASE_DIR, sim_name, f"{sim_name.replace('md_', '')}.dcd")
    if os.path.exists(dcd_file):
        size_mb = os.path.getsize(dcd_file) / (1024 * 1024)
        return size_mb
    return 0

def check_process():
    """실행 중인 OpenMM 프로세스 확인"""
    try:
        result = subprocess.run(['ps', 'aux'], capture_output=True, text=True)
        for line in result.stdout.split('\n'):
            if 'openmm_run.py' in line and 'tripod' in line.lower():
                return True
    except:
        pass
    return False

def get_log_tail(sim_name, lines=5):
    """로그 파일 마지막 몇 줄 확인"""
    log_file = os.path.join(BASE_DIR, sim_name, "run.log")
    if os.path.exists(log_file):
        try:
            with open(log_file, 'r') as f:
                all_lines = f.readlines()
                return ''.join(all_lines[-lines:])
        except:
            pass
    return "로그 없음"

def main():
    print("=" * 80)
    print("Tripod MD 시뮬레이션 진행상황")
    print("=" * 80)
    print()
    
    # 프로세스 확인
    is_running = check_process()
    if is_running:
        print("🟢 OpenMM 프로세스 실행 중")
    else:
        print("⚪ OpenMM 프로세스 없음 (완료 또는 대기 중)")
    print()
    
    # 각 시뮬레이션 상태 확인
    for sim_name in SIMULATIONS:
        print(f"📊 {sim_name}")
        print("-" * 80)
        
        dcd_size = get_dcd_size(sim_name)
        expected = EXPECTED_FRAMES[sim_name]
        
        if dcd_size > 0:
            print(f"  DCD 파일: {dcd_size:.1f} MB")
            
            # 대략적인 진행률 추정 (1 frame ≈ 0.2 MB for tripod)
            estimated_frames = int(dcd_size / 0.2)
            progress = min(100, (estimated_frames / expected) * 100)
            print(f"  예상 진행률: {progress:.1f}% ({estimated_frames}/{expected} frames)")
            
            # 진행 바
            bar_length = 50
            filled = int(bar_length * progress / 100)
            bar = '█' * filled + '░' * (bar_length - filled)
            print(f"  [{bar}] {progress:.1f}%")
        else:
            print(f"  상태: 시작 전 또는 초기화 중")
        
        print()
    
    print("=" * 80)
    print("💡 실시간 로그 확인:")
    print("   tail -f results/md_tripod_1ns/run.log")
    print("   tail -f results/md_tripod_10ns/run.log")
    print("   tail -f results/md_tripod_100ns/run.log")
    print("=" * 80)

if __name__ == "__main__":
    main()
