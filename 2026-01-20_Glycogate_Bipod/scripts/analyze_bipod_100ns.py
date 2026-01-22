#!/usr/bin/env python3
"""
Bipod 100ns MD 시뮬레이션 결과 분석
"""
import os
import numpy as np
from MDAnalysis import Universe
from MDAnalysis.lib.distances import calc_bonds

# 파일 경로
PSF = "/home/pjho3tr/projects/Drug/2026-01-20_Glycogate_Bipod/data/solution builder/openmm/step3_input.psf"
DCD = "/home/pjho3tr/projects/Drug/2026-01-20_Glycogate_Bipod/results/md_bipod_100ns/bipod_100ns.dcd"
OUTDIR = "/home/pjho3tr/projects/Drug/2026-01-20_Glycogate_Bipod/results/md_bipod_100ns"
ARM_CSV = os.path.join(OUTDIR, "bipod_arm_analysis.csv")
GLUCOSE_CSV = os.path.join(OUTDIR, "glucose_distance.csv")

SEL = "resname LIG and not name H*"
THR40 = 40.0
THR50 = 50.0

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    u = Universe(PSF, DCD)
    lig = u.select_atoms(SEL)
    
    print(f"Total LIG heavy atoms: {lig.n_atoms}")
    print(f"Total frames: {len(u.trajectory)}")

    # frame0에서 core와 두 개의 tip 선택
    u.trajectory[0]
    com = lig.center_of_mass()
    d_com = np.linalg.norm(lig.positions - com, axis=1)
    core_local = int(np.argmin(d_com))
    core = lig.atoms[core_local]

    d_core = np.linalg.norm(lig.positions - core.position, axis=1)
    sorted_indices = np.argsort(d_core)[::-1]
    tip1_local = sorted_indices[0]
    tip2_local = sorted_indices[1]
    
    tip1 = lig.atoms[tip1_local]
    tip2 = lig.atoms[tip2_local]

    print("\n=== Selected atoms ===")
    print(f"CORE: ix={core.ix} name={core.name}")
    print(f"TIP1: ix={tip1.ix} name={tip1.name}")
    print(f"TIP2: ix={tip2.ix} name={tip2.name}")

    core_ag = u.atoms[core.ix:core.ix+1]
    tip1_ag = u.atoms[tip1.ix:tip1.ix+1]
    tip2_ag = u.atoms[tip2.ix:tip2.ix+1]

    # Glucose 그룹 식별 (10ns와 동일한 방법)
    far_atoms_mask = d_core >= 80.0
    far_atoms_indices = np.where(far_atoms_mask)[0]
    if len(far_atoms_indices) < 10:
        far_atoms_indices = sorted_indices[:50]
    
    far_atoms = lig.atoms[far_atoms_indices]
    tip1_pos = tip1.position
    tip2_pos = tip2.position
    
    d_to_tip1 = np.linalg.norm(far_atoms.positions - tip1_pos, axis=1)
    d_to_tip2 = np.linalg.norm(far_atoms.positions - tip2_pos, axis=1)
    
    glucose1_mask = d_to_tip1 < d_to_tip2
    glucose2_mask = ~glucose1_mask
    
    glucose1_indices = far_atoms_indices[np.where(glucose1_mask)[0]]
    glucose2_indices = far_atoms_indices[np.where(glucose2_mask)[0]]
    
    glucose1 = lig.atoms[glucose1_indices]
    glucose2 = lig.atoms[glucose2_indices]
    
    print(f"\nGlucose 1: {len(glucose1)} atoms")
    print(f"Glucose 2: {len(glucose2)} atoms")

    # Trajectory 분석
    arm_rows = []
    glucose_rows = []
    arm1_vals = []
    arm2_vals = []
    max_vals = []
    rg_vals = []
    glucose_vals = []

    print("\nAnalyzing trajectory...")
    for i, ts in enumerate(u.trajectory):
        arm1 = float(calc_bonds(core_ag.positions, tip1_ag.positions, box=ts.dimensions)[0])
        arm2 = float(calc_bonds(core_ag.positions, tip2_ag.positions, box=ts.dimensions)[0])
        max_arm = max(arm1, arm2)
        rg = float(lig.radius_of_gyration())
        
        com1 = glucose1.center_of_mass()
        com2 = glucose2.center_of_mass()
        glucose_dist = np.linalg.norm(com1 - com2)
        
        t_ns = ts.time / 1000.0
        
        arm_rows.append((t_ns, arm1, arm2, max_arm, rg))
        glucose_rows.append((t_ns, glucose_dist))
        
        arm1_vals.append(arm1)
        arm2_vals.append(arm2)
        max_vals.append(max_arm)
        rg_vals.append(rg)
        glucose_vals.append(glucose_dist)
        
        if (i + 1) % 100 == 0:
            print(f"  Frame {i+1}/{len(u.trajectory)}: t={t_ns:.2f} ns")

    # CSV 저장
    with open(ARM_CSV, "w") as f:
        f.write("time_ns,arm1_A,arm2_A,max_arm_A,Rg_A\n")
        for row in arm_rows:
            f.write(f"{row[0]:.6f},{row[1]:.6f},{row[2]:.6f},{row[3]:.6f},{row[4]:.6f}\n")

    with open(GLUCOSE_CSV, "w") as f:
        f.write("time_ns,glucose_distance_A\n")
        for row in glucose_rows:
            f.write(f"{row[0]:.6f},{row[1]:.6f}\n")

    # 통계
    arm1_vals = np.array(arm1_vals)
    arm2_vals = np.array(arm2_vals)
    max_vals = np.array(max_vals)
    rg_vals = np.array(rg_vals)
    glucose_vals = np.array(glucose_vals)

    print("\n" + "=" * 80)
    print("RESULTS (100 ns)")
    print("=" * 80)
    print(f"Frames: {len(arm_rows)}")
    print()
    print("ARM 1 (Å):")
    print(f"  min/mean/max = {arm1_vals.min():.3f} / {arm1_vals.mean():.3f} / {arm1_vals.max():.3f}")
    print(f"  std = {arm1_vals.std():.3f}")
    print(f"  P(>= 40Å) = {(arm1_vals >= THR40).mean() * 100:.2f}%")
    print(f"  P(>= 50Å) = {(arm1_vals >= THR50).mean() * 100:.2f}%")
    print()
    print("ARM 2 (Å):")
    print(f"  min/mean/max = {arm2_vals.min():.3f} / {arm2_vals.mean():.3f} / {arm2_vals.max():.3f}")
    print(f"  std = {arm2_vals.std():.3f}")
    print(f"  P(>= 40Å) = {(arm2_vals >= THR40).mean() * 100:.2f}%")
    print(f"  P(>= 50Å) = {(arm2_vals >= THR50).mean() * 100:.2f}%")
    print()
    print("MAX ARM (Å):")
    print(f"  min/mean/max = {max_vals.min():.3f} / {max_vals.mean():.3f} / {max_vals.max():.3f}")
    print(f"  std = {max_vals.std():.3f}")
    print(f"  P(>= 40Å) = {(max_vals >= THR40).mean() * 100:.2f}%")
    print(f"  P(>= 50Å) = {(max_vals >= THR50).mean() * 100:.2f}%")
    print()
    print("Glucose Distance (Å):")
    print(f"  min/mean/max = {glucose_vals.min():.3f} / {glucose_vals.mean():.3f} / {glucose_vals.max():.3f}")
    print(f"  std = {glucose_vals.std():.3f}")
    print()
    print("Radius of Gyration (Å):")
    print(f"  min/mean/max = {rg_vals.min():.3f} / {rg_vals.mean():.3f} / {rg_vals.max():.3f}")
    print(f"  std = {rg_vals.std():.3f}")
    print("=" * 80)
    
    print(f"\n✅ Arm analysis: {ARM_CSV}")
    print(f"✅ Glucose distance: {GLUCOSE_CSV}")

if __name__ == "__main__":
    main()
