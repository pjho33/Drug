#!/usr/bin/env python
"""
Phase 2 CUDA 재실행 스크립트
Tripod, Glucose, BNG 3개 리간드 100ns 시뮬레이션 (GPU 가속)
"""
import os
import sys
import subprocess
import shutil
from datetime import datetime

def log(msg):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}")

def run_simulation(receptor, ligand_sdf, out_prefix, total_ns=100.0):
    """Run a single simulation using 03_run_production.py with CUDA"""
    script = "/home/pjho3/projects/Drug/scripts/03_run_production.py"
    cmd = [
        sys.executable, script,
        receptor, ligand_sdf, out_prefix,
        "--total_ns", str(total_ns),
        "--platform", "CUDA",
        "--device_index", "0",
        "--precision", "mixed"
    ]
    log(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=os.path.dirname(out_prefix))
    return result.returncode == 0

def main():
    base_dir = "/home/pjho3/projects/Drug/results"
    rep1_dir = os.path.join(base_dir, "phase2_rep1")
    
    # 기존 결과 백업
    backup_dir = os.path.join(base_dir, "phase2_rep1_backup_cpu")
    if not os.path.exists(backup_dir):
        log(f"Backing up old results to {backup_dir}")
        os.makedirs(backup_dir, exist_ok=True)
        for f in os.listdir(rep1_dir):
            if f.startswith("prod_") and not f.endswith(":Zone.Identifier"):
                src = os.path.join(rep1_dir, f)
                dst = os.path.join(backup_dir, f)
                if os.path.isfile(src):
                    shutil.move(src, dst)
                    log(f"  Moved: {f}")
    
    receptor = os.path.join(rep1_dir, "cleaned_receptor.pdb")
    ligands = {
        "tripod": os.path.join(rep1_dir, "dock_tripod", "complex_0", "rank1.sdf"),
        "glucose": os.path.join(rep1_dir, "dock_glucose", "complex_0", "rank1.sdf"),
        "bng": os.path.join(rep1_dir, "dock_bng", "complex_0", "rank1.sdf"),
    }
    
    log("===== Phase 2 CUDA Re-run (RTX 3090) =====")
    log("Expected time: ~30min - 2hr per simulation")
    
    for name in ["tripod", "glucose", "bng"]:
        out_prefix = os.path.join(rep1_dir, f"prod_{name}_rep1")
        final_pdb = out_prefix + "_final.pdb"
        
        # 이미 완료된 경우 스킵
        if os.path.exists(final_pdb):
            log(f"Skipping {name} - already completed")
            continue
            
        log(f"Starting {name} (100ns)...")
        success = run_simulation(receptor, ligands[name], out_prefix)
        if success:
            log(f"{name} completed!")
        else:
            log(f"{name} FAILED!")
            return 1
    
    log("===== All Phase 2 rep1 simulations completed! =====")
    return 0

if __name__ == "__main__":
    sys.exit(main())
