#!/usr/bin/env python3
"""
Tripod 100ns MD 시뮬레이션 결과 분석
- 3개 arm 길이 측정
- Glucose 간 거리 측정
- Radius of gyration
"""
import os
import numpy as np
from MDAnalysis import Universe
from MDAnalysis.lib.distances import calc_bonds

# 파일 경로
PSF = "/home/pjho3tr/projects/Drug/2026-01-22_Glycogate_Tripod/data/solution_builder/openmm/step3_input.psf"
DCD = "/home/pjho3tr/projects/Drug/2026-01-22_Glycogate_Tripod/results/md_tripod_100ns/tripod_100ns.dcd"
OUTDIR = "/home/pjho3tr/projects/Drug/2026-01-22_Glycogate_Tripod/results/md_tripod_100ns"

SEL = "resname LIG and not name H*"

def main():
    u = Universe(PSF, DCD)
    lig = u.select_atoms(SEL)
    
    print(f"Total LIG heavy atoms: {lig.n_atoms}")
    print(f"Total frames: {len(u.trajectory)}")
    
    # Frame 0에서 core와 3개 tip 선택
    u.trajectory[0]
    com = lig.center_of_mass()
    d_com = np.linalg.norm(lig.positions - com, axis=1)
    core_local = int(np.argmin(d_com))
    core = lig.atoms[core_local]
    
    # 가장 먼 3개 원자를 tip으로 선택
    d_core = np.linalg.norm(lig.positions - core.position, axis=1)
    sorted_indices = np.argsort(d_core)[::-1]
    
    tip1_local = sorted_indices[0]
    tip2_local = sorted_indices[1]
    tip3_local = sorted_indices[2]
    
    tip1 = lig.atoms[tip1_local]
    tip2 = lig.atoms[tip2_local]
    tip3 = lig.atoms[tip3_local]
    
    print("\n=== Selected atoms ===")
    print(f"CORE: ix={core.ix} name={core.name}")
    print(f"TIP1: ix={tip1.ix} name={tip1.name}")
    print(f"TIP2: ix={tip2.ix} name={tip2.name}")
    print(f"TIP3: ix={tip3.ix} name={tip3.name}")
    
    core_ag = u.atoms[core.ix:core.ix+1]
    tip1_ag = u.atoms[tip1.ix:tip1.ix+1]
    tip2_ag = u.atoms[tip2.ix:tip2.ix+1]
    tip3_ag = u.atoms[tip3.ix:tip3.ix+1]
    
    # Glucose 그룹 식별 (각 tip 주변)
    far_atoms_mask = d_core >= 80.0
    far_atoms_indices = np.where(far_atoms_mask)[0]
    if len(far_atoms_indices) < 15:
        far_atoms_indices = sorted_indices[:100]
    
    far_atoms = lig.atoms[far_atoms_indices]
    tip1_pos = tip1.position
    tip2_pos = tip2.position
    tip3_pos = tip3.position
    
    d_to_tip1 = np.linalg.norm(far_atoms.positions - tip1_pos, axis=1)
    d_to_tip2 = np.linalg.norm(far_atoms.positions - tip2_pos, axis=1)
    d_to_tip3 = np.linalg.norm(far_atoms.positions - tip3_pos, axis=1)
    
    # 각 원자를 가장 가까운 tip에 할당
    min_distances = np.minimum(np.minimum(d_to_tip1, d_to_tip2), d_to_tip3)
    glucose1_mask = (d_to_tip1 == min_distances)
    glucose2_mask = (d_to_tip2 == min_distances)
    glucose3_mask = (d_to_tip3 == min_distances)
    
    glucose1_indices = far_atoms_indices[np.where(glucose1_mask)[0]]
    glucose2_indices = far_atoms_indices[np.where(glucose2_mask)[0]]
    glucose3_indices = far_atoms_indices[np.where(glucose3_mask)[0]]
    
    glucose1 = lig.atoms[glucose1_indices]
    glucose2 = lig.atoms[glucose2_indices]
    glucose3 = lig.atoms[glucose3_indices]
    
    print(f"\nGlucose 1: {len(glucose1)} atoms")
    print(f"Glucose 2: {len(glucose2)} atoms")
    print(f"Glucose 3: {len(glucose3)} atoms")
    
    # Trajectory 분석
    arm_rows = []
    glucose_rows = []
    
    arm1_vals = []
    arm2_vals = []
    arm3_vals = []
    max_vals = []
    rg_vals = []
    
    glucose12_vals = []
    glucose13_vals = []
    glucose23_vals = []
    
    print("\nAnalyzing trajectory...")
    for i, ts in enumerate(u.trajectory):
        arm1 = float(calc_bonds(core_ag.positions, tip1_ag.positions, box=ts.dimensions)[0])
        arm2 = float(calc_bonds(core_ag.positions, tip2_ag.positions, box=ts.dimensions)[0])
        arm3 = float(calc_bonds(core_ag.positions, tip3_ag.positions, box=ts.dimensions)[0])
        max_arm = max(arm1, arm2, arm3)
        rg = float(lig.radius_of_gyration())
        
        com1 = glucose1.center_of_mass()
        com2 = glucose2.center_of_mass()
        com3 = glucose3.center_of_mass()
        
        glucose12 = np.linalg.norm(com1 - com2)
        glucose13 = np.linalg.norm(com1 - com3)
        glucose23 = np.linalg.norm(com2 - com3)
        
        t_ns = ts.time / 1000.0
        
        arm_rows.append((t_ns, arm1, arm2, arm3, max_arm, rg))
        glucose_rows.append((t_ns, glucose12, glucose13, glucose23))
        
        arm1_vals.append(arm1)
        arm2_vals.append(arm2)
        arm3_vals.append(arm3)
        max_vals.append(max_arm)
        rg_vals.append(rg)
        
        glucose12_vals.append(glucose12)
        glucose13_vals.append(glucose13)
        glucose23_vals.append(glucose23)
        
        if (i + 1) % 100 == 0:
            print(f"  Frame {i+1}/{len(u.trajectory)}: t={t_ns:.3f} ns")
    
    # CSV 저장
    arm_csv = os.path.join(OUTDIR, "tripod_arm_analysis.csv")
    with open(arm_csv, "w") as f:
        f.write("time_ns,arm1_A,arm2_A,arm3_A,max_arm_A,Rg_A\n")
        for row in arm_rows:
            f.write(f"{row[0]:.6f},{row[1]:.6f},{row[2]:.6f},{row[3]:.6f},{row[4]:.6f},{row[5]:.6f}\n")
    
    glucose_csv = os.path.join(OUTDIR, "glucose_distances.csv")
    with open(glucose_csv, "w") as f:
        f.write("time_ns,glucose12_A,glucose13_A,glucose23_A\n")
        for row in glucose_rows:
            f.write(f"{row[0]:.6f},{row[1]:.6f},{row[2]:.6f},{row[3]:.6f}\n")
    
    # 통계
    arm1_vals = np.array(arm1_vals)
    arm2_vals = np.array(arm2_vals)
    arm3_vals = np.array(arm3_vals)
    max_vals = np.array(max_vals)
    rg_vals = np.array(rg_vals)
    
    glucose12_vals = np.array(glucose12_vals)
    glucose13_vals = np.array(glucose13_vals)
    glucose23_vals = np.array(glucose23_vals)
    
    print("\n" + "=" * 80)
    print("RESULTS (100 ns)")
    print("=" * 80)
    print(f"Frames: {len(arm_rows)}")
    print()
    print("ARM 1 (Å):")
    print(f"  mean ± std = {arm1_vals.mean():.2f} ± {arm1_vals.std():.2f}")
    print(f"  min / max = {arm1_vals.min():.2f} / {arm1_vals.max():.2f}")
    print()
    print("ARM 2 (Å):")
    print(f"  mean ± std = {arm2_vals.mean():.2f} ± {arm2_vals.std():.2f}")
    print(f"  min / max = {arm2_vals.min():.2f} / {arm2_vals.max():.2f}")
    print()
    print("ARM 3 (Å):")
    print(f"  mean ± std = {arm3_vals.mean():.2f} ± {arm3_vals.std():.2f}")
    print(f"  min / max = {arm3_vals.min():.2f} / {arm3_vals.max():.2f}")
    print()
    print("MAX ARM (Å):")
    print(f"  mean ± std = {max_vals.mean():.2f} ± {max_vals.std():.2f}")
    print(f"  min / max = {max_vals.min():.2f} / {max_vals.max():.2f}")
    print(f"  P(>= 40Å) = {(max_vals >= 40).mean() * 100:.1f}%")
    print(f"  P(>= 50Å) = {(max_vals >= 50).mean() * 100:.1f}%")
    print()
    print("Glucose Distances (Å):")
    print(f"  G1-G2: {glucose12_vals.mean():.2f} ± {glucose12_vals.std():.2f}")
    print(f"  G1-G3: {glucose13_vals.mean():.2f} ± {glucose13_vals.std():.2f}")
    print(f"  G2-G3: {glucose23_vals.mean():.2f} ± {glucose23_vals.std():.2f}")
    print()
    print("Radius of Gyration (Å):")
    print(f"  mean ± std = {rg_vals.mean():.2f} ± {rg_vals.std():.2f}")
    print(f"  min / max = {rg_vals.min():.2f} / {rg_vals.max():.2f}")
    print("=" * 80)
    
    print(f"\n✅ Arm analysis: {arm_csv}")
    print(f"✅ Glucose distances: {glucose_csv}")

if __name__ == "__main__":
    main()
