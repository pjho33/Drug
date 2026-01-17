#!/usr/bin/env python
"""
Phase 2 자동 실행 스크립트
rep2 (glucose, bng) + rep3 (tripod, glucose, bng) 순차 실행
"""
import os
import sys
import subprocess
from datetime import datetime

def log(msg):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}")

def run_simulation(receptor, ligand_sdf, out_prefix, total_ns=100.0):
    """Run a single simulation using 03_run_production.py"""
    script = "/home/pjho3/projects/Drug/scripts/03_run_production.py"
    cmd = [
        sys.executable, script,
        receptor, ligand_sdf, out_prefix,
        "--total_ns", str(total_ns)
    ]
    log(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=os.path.dirname(out_prefix))
    return result.returncode == 0

def main():
    base_dir = "/home/pjho3/projects/Drug/results"
    rep1_dir = os.path.join(base_dir, "phase2_rep1")
    
    receptor = os.path.join(rep1_dir, "cleaned_receptor.pdb")
    ligands = {
        "tripod": os.path.join(rep1_dir, "dock_tripod", "complex_0", "rank1.sdf"),
        "glucose": os.path.join(rep1_dir, "dock_glucose", "complex_0", "rank1.sdf"),
        "bng": os.path.join(rep1_dir, "dock_bng", "complex_0", "rank1.sdf"),
    }
    
    log("===== Phase 2 Auto Run Script =====")
    
    # Rep2 - glucose and bng (tripod is already running)
    rep2_dir = os.path.join(base_dir, "phase2_rep2")
    os.makedirs(rep2_dir, exist_ok=True)
    
    for name in ["glucose", "bng"]:
        out_prefix = os.path.join(rep2_dir, f"prod_{name}_rep2")
        final_pdb = out_prefix + "_final.pdb"
        if os.path.exists(final_pdb):
            log(f"Skipping {name} rep2 - already completed")
            continue
        log(f"Starting {name} rep2...")
        success = run_simulation(receptor, ligands[name], out_prefix)
        if success:
            log(f"{name} rep2 completed!")
        else:
            log(f"{name} rep2 FAILED!")
    
    # Rep3 - all three
    rep3_dir = os.path.join(base_dir, "phase2_rep3")
    os.makedirs(rep3_dir, exist_ok=True)
    
    for name in ["tripod", "glucose", "bng"]:
        out_prefix = os.path.join(rep3_dir, f"prod_{name}_rep3")
        final_pdb = out_prefix + "_final.pdb"
        if os.path.exists(final_pdb):
            log(f"Skipping {name} rep3 - already completed")
            continue
        log(f"Starting {name} rep3...")
        success = run_simulation(receptor, ligands[name], out_prefix)
        if success:
            log(f"{name} rep3 completed!")
        else:
            log(f"{name} rep3 FAILED!")
    
    log("===== All Phase 2 simulations completed! =====")

if __name__ == "__main__":
    main()
