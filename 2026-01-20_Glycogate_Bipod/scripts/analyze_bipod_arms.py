#!/usr/bin/env python3
"""
Bipod 리간드 분석: 두 개의 arm 길이 측정

- core: frame0에서 LIG heavy-atom COM에 가장 가까운 원자
- tip1, tip2: frame0에서 core로부터 가장 먼 heavy atom 2개 (두 팔 끝)
- 출력: CSV (time_ns, arm1_A, arm2_A, max_arm_A, Rg_A)
- 4nm/5nm(40/50Å) 이상 비율 출력
"""
import os
import numpy as np
from MDAnalysis import Universe
from MDAnalysis.lib.distances import calc_bonds

# 파일 경로
PSF = "/home/pjho3tr/projects/Drug/2026-01-20_Glycogate_Bipod/data/solution builder/openmm/step3_input.psf"
DCD = "/home/pjho3tr/projects/Drug/2026-01-20_Glycogate_Bipod/results/md_bipod_10ns/bipod_10ns.dcd"
OUTDIR = "/home/pjho3tr/projects/Drug/2026-01-20_Glycogate_Bipod/results/md_bipod_10ns"
OUTCSV = os.path.join(OUTDIR, "bipod_arm_analysis.csv")

SEL = "resname LIG and not name H*"
THR40 = 40.0
THR50 = 50.0

def info(a):
    return f"ix={a.ix} name={a.name} type={a.type} resid={a.resid} resname={a.resname}"

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    u = Universe(PSF, DCD)
    lig = u.select_atoms(SEL)
    if lig.n_atoms == 0:
        raise SystemExit(f"ERROR: empty selection: {SEL}")

    print(f"Total LIG heavy atoms: {lig.n_atoms}")

    # frame0에서 core와 두 개의 tip 선택
    u.trajectory[0]
    com = lig.center_of_mass()
    d_com = np.linalg.norm(lig.positions - com, axis=1)
    core_local = int(np.argmin(d_com))
    core = lig.atoms[core_local]

    # core로부터 거리 계산
    d_core = np.linalg.norm(lig.positions - core.position, axis=1)
    
    # 가장 먼 두 개의 원자를 tip으로 선택
    sorted_indices = np.argsort(d_core)[::-1]  # 내림차순
    tip1_local = sorted_indices[0]
    tip2_local = sorted_indices[1]
    
    tip1 = lig.atoms[tip1_local]
    tip2 = lig.atoms[tip2_local]

    print("\n=== Selected atoms ===")
    print("CORE:", info(core))
    print("TIP1:", info(tip1), f"distance={d_core[tip1_local]:.2f} Å")
    print("TIP2:", info(tip2), f"distance={d_core[tip2_local]:.2f} Å")

    core_ag = u.atoms[core.ix:core.ix+1]
    tip1_ag = u.atoms[tip1.ix:tip1.ix+1]
    tip2_ag = u.atoms[tip2.ix:tip2.ix+1]

    # trajectory 분석
    rows = []
    arm1_vals = []
    arm2_vals = []
    max_vals = []
    rg_vals = []

    print("\nAnalyzing trajectory...")
    for i, ts in enumerate(u.trajectory):
        arm1 = float(calc_bonds(core_ag.positions, tip1_ag.positions, box=ts.dimensions)[0])
        arm2 = float(calc_bonds(core_ag.positions, tip2_ag.positions, box=ts.dimensions)[0])
        max_arm = max(arm1, arm2)
        rg = float(lig.radius_of_gyration())
        t_ns = ts.time / 1000.0
        
        rows.append((t_ns, arm1, arm2, max_arm, rg))
        arm1_vals.append(arm1)
        arm2_vals.append(arm2)
        max_vals.append(max_arm)
        rg_vals.append(rg)
        
        if (i + 1) % 10 == 0:
            print(f"  Frame {i+1}/{len(u.trajectory)}: t={t_ns:.2f} ns, arm1={arm1:.2f} Å, arm2={arm2:.2f} Å")

    # CSV 저장
    with open(OUTCSV, "w") as f:
        f.write("time_ns,arm1_A,arm2_A,max_arm_A,Rg_A\n")
        for row in rows:
            f.write(f"{row[0]:.6f},{row[1]:.6f},{row[2]:.6f},{row[3]:.6f},{row[4]:.6f}\n")

    # 통계 출력
    arm1_vals = np.array(arm1_vals)
    arm2_vals = np.array(arm2_vals)
    max_vals = np.array(max_vals)
    rg_vals = np.array(rg_vals)

    print("\n" + "=" * 80)
    print("RESULTS")
    print("=" * 80)
    print(f"Output: {OUTCSV}")
    print(f"Frames: {len(rows)}")
    print()
    print("ARM 1 (Å):")
    print(f"  min/mean/max = {arm1_vals.min():.3f} / {arm1_vals.mean():.3f} / {arm1_vals.max():.3f}")
    print(f"  P(>= 40Å / 4nm) = {(arm1_vals >= THR40).mean() * 100:.2f}%")
    print(f"  P(>= 50Å / 5nm) = {(arm1_vals >= THR50).mean() * 100:.2f}%")
    print()
    print("ARM 2 (Å):")
    print(f"  min/mean/max = {arm2_vals.min():.3f} / {arm2_vals.mean():.3f} / {arm2_vals.max():.3f}")
    print(f"  P(>= 40Å / 4nm) = {(arm2_vals >= THR40).mean() * 100:.2f}%")
    print(f"  P(>= 50Å / 5nm) = {(arm2_vals >= THR50).mean() * 100:.2f}%")
    print()
    print("MAX ARM (Å):")
    print(f"  min/mean/max = {max_vals.min():.3f} / {max_vals.mean():.3f} / {max_vals.max():.3f}")
    print(f"  P(>= 40Å / 4nm) = {(max_vals >= THR40).mean() * 100:.2f}%")
    print(f"  P(>= 50Å / 5nm) = {(max_vals >= THR50).mean() * 100:.2f}%")
    print()
    print("Radius of Gyration (Å):")
    print(f"  min/mean/max = {rg_vals.min():.3f} / {rg_vals.mean():.3f} / {rg_vals.max():.3f}")
    print("=" * 80)

if __name__ == "__main__":
    main()
