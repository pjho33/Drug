#!/usr/bin/env python3
"""
Bipod 리간드: 두 L-glucose 사이의 거리 측정 (개선 버전)

두 팔 끝의 glucose를 클러스터링으로 식별하여 거리 측정
"""
import os
import numpy as np
from MDAnalysis import Universe
from scipy.spatial.distance import cdist

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

    # frame0에서 두 glucose 그룹 식별
    u.trajectory[0]
    com = lig.center_of_mass()
    d_com = np.linalg.norm(lig.positions - com, axis=1)
    core_local = int(np.argmin(d_com))
    core = lig.atoms[core_local]

    # core로부터 거리 계산
    d_core = np.linalg.norm(lig.positions - core.position, axis=1)
    
    # 가장 먼 두 개의 원자를 각 glucose의 대표 원자로 선택
    sorted_indices = np.argsort(d_core)[::-1]
    tip1_local = sorted_indices[0]
    tip2_local = sorted_indices[1]
    
    tip1 = lig.atoms[tip1_local]
    tip2 = lig.atoms[tip2_local]
    
    tip1_pos = tip1.position
    tip2_pos = tip2.position

    print("\n=== Tip atoms (glucose representatives) ===")
    print(f"TIP1: ix={tip1.ix} name={tip1.name} type={tip1.type} pos={tip1_pos}")
    print(f"TIP2: ix={tip2.ix} name={tip2.name} type={tip2.type} pos={tip2_pos}")
    print(f"Distance between tips: {np.linalg.norm(tip1_pos - tip2_pos):.2f} Å")

    # 각 원자를 tip1 또는 tip2에 가까운 쪽으로 분류
    # core에서 멀리 떨어진 원자들만 고려 (core로부터 80 Å 이상)
    far_atoms_mask = d_core >= 80.0
    far_atoms_indices = np.where(far_atoms_mask)[0]
    far_atoms = lig.atoms[far_atoms_indices]
    
    print(f"\nAtoms far from core (>= 80 Å): {len(far_atoms)}")
    
    if len(far_atoms) < 10:
        print("WARNING: Too few atoms far from core. Adjusting threshold...")
        # 상위 50개 원자 선택
        far_atoms_indices = sorted_indices[:50]
        far_atoms = lig.atoms[far_atoms_indices]
        print(f"Using top 50 farthest atoms instead: {len(far_atoms)}")
    
    # 각 원자를 tip1 또는 tip2에 할당
    d_to_tip1 = np.linalg.norm(far_atoms.positions - tip1_pos, axis=1)
    d_to_tip2 = np.linalg.norm(far_atoms.positions - tip2_pos, axis=1)
    
    glucose1_mask = d_to_tip1 < d_to_tip2
    glucose2_mask = ~glucose1_mask
    
    glucose1_local_indices = np.where(glucose1_mask)[0]
    glucose2_local_indices = np.where(glucose2_mask)[0]
    
    glucose1_indices = far_atoms_indices[glucose1_local_indices]
    glucose2_indices = far_atoms_indices[glucose2_local_indices]
    
    glucose1 = lig.atoms[glucose1_indices]
    glucose2 = lig.atoms[glucose2_indices]
    
    print(f"\nGlucose 1 group: {len(glucose1)} atoms")
    print(f"Glucose 2 group: {len(glucose2)} atoms")
    
    if len(glucose1) == 0 or len(glucose2) == 0:
        raise SystemExit("ERROR: Could not identify glucose groups")
    
    # 초기 프레임에서 두 glucose COM 사이 거리 확인
    com1_init = glucose1.center_of_mass()
    com2_init = glucose2.center_of_mass()
    init_distance = np.linalg.norm(com1_init - com2_init)
    print(f"\nInitial distance between glucose COMs: {init_distance:.2f} Å")

    # trajectory 분석
    rows = []
    distance_vals = []

    print("\nAnalyzing trajectory...")
    for i, ts in enumerate(u.trajectory):
        # 각 glucose의 COM 계산
        com1 = glucose1.center_of_mass()
        com2 = glucose2.center_of_mass()
        
        # 두 COM 사이의 거리
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
