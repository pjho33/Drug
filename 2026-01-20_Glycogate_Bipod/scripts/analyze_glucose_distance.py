#!/usr/bin/env python3
"""
Bipod 리간드: 두 L-glucose 사이의 거리 측정

- 각 팔 끝의 glucose를 식별
- 두 glucose의 COM(center of mass) 사이 거리 측정
- 출력: CSV (time_ns, glucose_distance_A)
"""
import os
import numpy as np
from MDAnalysis import Universe
from MDAnalysis.lib.distances import calc_bonds

# 파일 경로
PSF = "/home/pjho3tr/projects/Drug/2026-01-20_Glycogate_Bipod/data/solution builder/openmm/step3_input.psf"
DCD = "/home/pjho3tr/projects/Drug/2026-01-20_Glycogate_Bipod/results/md_bipod_10ns/bipod_10ns.dcd"
OUTDIR = "/home/pjho3tr/projects/Drug/2026-01-20_Glycogate_Bipod/results/md_bipod_10ns"
OUTCSV = os.path.join(OUTDIR, "glucose_distance.csv")

SEL = "resname LIG and not name H*"

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    u = Universe(PSF, DCD)
    lig = u.select_atoms(SEL)
    if lig.n_atoms == 0:
        raise SystemExit(f"ERROR: empty selection: {SEL}")

    print(f"Total LIG heavy atoms: {lig.n_atoms}")

    # frame0에서 core와 두 개의 tip 선택 (이전 분석과 동일한 방법)
    u.trajectory[0]
    com = lig.center_of_mass()
    d_com = np.linalg.norm(lig.positions - com, axis=1)
    core_local = int(np.argmin(d_com))
    core = lig.atoms[core_local]

    # core로부터 거리 계산
    d_core = np.linalg.norm(lig.positions - core.position, axis=1)
    
    # 가장 먼 두 개의 원자를 tip으로 선택
    sorted_indices = np.argsort(d_core)[::-1]
    tip1_local = sorted_indices[0]
    tip2_local = sorted_indices[1]
    
    tip1 = lig.atoms[tip1_local]
    tip2 = lig.atoms[tip2_local]

    print("\n=== Tip atoms (glucose ends) ===")
    print(f"TIP1: ix={tip1.ix} name={tip1.name} resid={tip1.resid}")
    print(f"TIP2: ix={tip2.ix} name={tip2.name} resid={tip2.resid}")

    # 각 tip 주변의 glucose 원자들을 선택 (tip으로부터 10 Å 이내)
    # glucose는 보통 6-7개의 heavy atom으로 구성
    tip1_pos = tip1.position
    tip2_pos = tip2.position
    
    # tip 주변 원자 찾기 (glucose 구조)
    d_tip1 = np.linalg.norm(lig.positions - tip1_pos, axis=1)
    d_tip2 = np.linalg.norm(lig.positions - tip2_pos, axis=1)
    
    # 각 tip으로부터 10 Å 이내의 원자들을 glucose로 간주
    glucose1_mask = d_tip1 <= 10.0
    glucose2_mask = d_tip2 <= 10.0
    
    glucose1_indices = np.where(glucose1_mask)[0]
    glucose2_indices = np.where(glucose2_mask)[0]
    
    # AtomGroup 생성
    glucose1 = lig.atoms[glucose1_indices]
    glucose2 = lig.atoms[glucose2_indices]
    
    print(f"\nGlucose 1: {len(glucose1)} atoms")
    print(f"Glucose 2: {len(glucose2)} atoms")
    
    if len(glucose1) == 0 or len(glucose2) == 0:
        raise SystemExit("ERROR: Could not identify glucose groups")

    # trajectory 분석
    rows = []
    distance_vals = []

    print("\nAnalyzing trajectory...")
    for i, ts in enumerate(u.trajectory):
        # 각 glucose의 COM 계산
        com1 = glucose1.center_of_mass()
        com2 = glucose2.center_of_mass()
        
        # 두 COM 사이의 거리 (PBC 고려)
        distance = np.linalg.norm(com1 - com2)
        
        t_ns = ts.time / 1000.0
        rows.append((t_ns, distance))
        distance_vals.append(distance)
        
        if (i + 1) % 10 == 0:
            print(f"  Frame {i+1}/{len(u.trajectory)}: t={t_ns:.2f} ns, distance={distance:.2f} Å")

    # CSV 저장
    with open(OUTCSV, "w") as f:
        f.write("time_ns,glucose_distance_A\n")
        for row in rows:
            f.write(f"{row[0]:.6f},{row[1]:.6f}\n")

    # 통계 출력
    distance_vals = np.array(distance_vals)

    print("\n" + "=" * 80)
    print("RESULTS: Distance between two L-glucose groups")
    print("=" * 80)
    print(f"Output: {OUTCSV}")
    print(f"Frames: {len(rows)}")
    print()
    print("Glucose-Glucose Distance (Å):")
    print(f"  min/mean/max = {distance_vals.min():.3f} / {distance_vals.mean():.3f} / {distance_vals.max():.3f}")
    print(f"  std = {distance_vals.std():.3f}")
    print()
    print("Distance distribution:")
    print(f"  < 50 Å  : {(distance_vals < 50).sum()} frames ({(distance_vals < 50).mean() * 100:.1f}%)")
    print(f"  50-100 Å: {((distance_vals >= 50) & (distance_vals < 100)).sum()} frames ({((distance_vals >= 50) & (distance_vals < 100)).mean() * 100:.1f}%)")
    print(f"  100-150 Å: {((distance_vals >= 100) & (distance_vals < 150)).sum()} frames ({((distance_vals >= 100) & (distance_vals < 150)).mean() * 100:.1f}%)")
    print(f"  >= 150 Å: {(distance_vals >= 150).sum()} frames ({(distance_vals >= 150).mean() * 100:.1f}%)")
    print("=" * 80)

if __name__ == "__main__":
    main()
