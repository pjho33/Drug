#!/usr/bin/env python3
"""
본 분석 스크립트: 한 다리 고정 조건에서 나머지 두 다리 도달거리/동시성 계산
Tripod Multivalency Analysis with Anchor Constraint
"""

import sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

# ==========================
# USER SETTINGS (EDIT HERE)
# ==========================

# 1) Target region (vestibule 주변 residue set)
#    예시: chain/segid에 맞게 수정하세요.
#    00_inspect_selections.py 결과를 참고하여 수정
TARGET_SEL = "protein and resid 50-80 150-190 280-320"  # <-- 반드시 수정 권장

# 2) Tripod 3 arms: 말단 glucose(또는 말단 그룹) selection 3개
#    가장 쉬운 방식: resid로 찍기
#    00_inspect_selections.py 결과에서 Glucose resid 3개를 찾아서 입력
ARM_SELS = [
    "resname GLC and resid 900",  # arm A
    "resname GLC and resid 901",  # arm B
    "resname GLC and resid 902",  # arm C
]

# 3) 각 arm에서 "말단 대표 원자" 선택
#    COM: sugar ring 전체 질량중심
#    ATOM: 특정 원자 하나 (예: "name O1")
ARM_REP_MODE = "COM"  # "COM" or "ATOM"
ARM_REP_ATOMSEL = "name O1"  # ARM_REP_MODE="ATOM"일 때 사용

# 4) Anchor/Reach cutoffs (Å)
ANCHOR_CUTOFF = 3.5  # arm이 타겟에 붙었다고 보는 기준
REACH_CUTOFFS = [3.5, 5.0, 8.0]  # 여러 임계값으로 스윕 추천

# 5) Output
OUT_NPZ = "tripod_anchor_reach.npz"

# ==========================


def rep_point(group, mode="COM", atomsel="name O1"):
    """Return 3D point of group: COM or a specific atom."""
    if len(group) == 0:
        raise ValueError("Empty atomgroup in rep_point()")
    if mode.upper() == "COM":
        return group.center_of_mass()
    elif mode.upper() == "ATOM":
        sub = group.select_atoms(atomsel)
        if len(sub) != 1:
            raise ValueError(f"ARM_REP_ATOMSEL should match exactly 1 atom, got {len(sub)}")
        return sub.positions[0]
    else:
        raise ValueError("mode must be COM or ATOM")


def min_dist_point_to_target(point_xyz, target_positions):
    """Min distance from a point to a set of target atoms positions (Å)."""
    # point_xyz shape (3,), target_positions shape (N,3)
    d = np.linalg.norm(target_positions - point_xyz[None, :], axis=1)
    return float(d.min())


def main(top, traj):
    print("\n" + "="*80)
    print("Tripod Anchor-Reach Analysis")
    print("="*80 + "\n")
    
    u = mda.Universe(top, traj)

    target = u.select_atoms(TARGET_SEL)
    if len(target) == 0:
        raise ValueError(f"TARGET_SEL matched 0 atoms: {TARGET_SEL}")

    arms = [u.select_atoms(s) for s in ARM_SELS]
    for i, ag in enumerate(arms):
        if len(ag) == 0:
            raise ValueError(f"ARM_SELS[{i}] matched 0 atoms: {ARM_SELS[i]}")

    n_frames = len(u.trajectory)
    print(f"Frames: {n_frames}")
    print(f"Target atoms: {len(target)}")
    for i, ag in enumerate(arms):
        print(f"Arm {i} atoms: {len(ag)} | sel: {ARM_SELS[i]}")
    print()

    # Store per-frame: per-arm min distance to target
    arm_to_target = np.zeros((n_frames, 3), dtype=np.float32)

    # Also store which arm is anchored (min distance smallest), and whether any anchored under cutoff
    anchor_arm = np.full(n_frames, -1, dtype=np.int32)
    anchor_dist = np.full(n_frames, np.nan, dtype=np.float32)

    # For each reach cutoff: store counts of (>=1 other arm reach) and (>=2 other arms reach) conditional on anchor
    reach_stats = {rc: {"anchor_frames": 0, "ge1": 0, "ge2": 0, "ge3": 0} for rc in REACH_CUTOFFS}

    # Dwell time: consecutive frames where >=2 arms reach under cutoff, within anchor frames
    ge2_series = {rc: np.zeros(n_frames, dtype=bool) for rc in REACH_CUTOFFS}
    anchor_series = np.zeros(n_frames, dtype=bool)

    print("Processing frames...")
    for fi, ts in enumerate(u.trajectory):
        if (fi + 1) % 1000 == 0:
            print(f"  Frame {fi+1}/{n_frames}")
            
        tgt_pos = target.positions  # (N,3)

        # per arm distance
        dists = []
        for ai, arm_ag in enumerate(arms):
            p = rep_point(arm_ag, mode=ARM_REP_MODE, atomsel=ARM_REP_ATOMSEL)
            d = min_dist_point_to_target(p, tgt_pos)
            arm_to_target[fi, ai] = d
            dists.append(d)

        dists = np.array(dists, dtype=float)
        best_arm = int(dists.argmin())
        best_d = float(dists[best_arm])

        anchor_arm[fi] = best_arm
        anchor_dist[fi] = best_d

        if best_d <= ANCHOR_CUTOFF:
            anchor_series[fi] = True

            for rc in REACH_CUTOFFS:
                reach_stats[rc]["anchor_frames"] += 1

                # Count how many arms are "reaching" under rc
                reaches = (dists <= rc)
                n_reach = int(reaches.sum())

                # Conditional: anchor is true already. Evaluate total multivalency.
                if n_reach >= 2:
                    reach_stats[rc]["ge1"] += 1
                if n_reach >= 3:
                    reach_stats[rc]["ge2"] += 1
                if n_reach >= 3:
                    reach_stats[rc]["ge3"] += 1

                ge2_series[rc][fi] = (n_reach >= 3)

    # Summaries
    print("\n" + "="*80)
    print("=== Anchor Summary ===")
    print("="*80)
    anchor_frac = anchor_series.mean()
    print(f"Anchor frames (best arm within {ANCHOR_CUTOFF}Å): {anchor_series.sum()} / {n_frames} = {anchor_frac:.3f}")

    print("\n" + "="*80)
    print("=== Conditional Reach Stats (given Anchor=True) ===")
    print("="*80)
    for rc in REACH_CUTOFFS:
        af = reach_stats[rc]["anchor_frames"]
        if af == 0:
            print(f"[rc={rc}Å] No anchor frames -> cannot compute conditional stats.")
            continue
        p_ge2arms = reach_stats[rc]["ge1"] / af     # at least 2 arms total reach
        p_ge3arms = reach_stats[rc]["ge2"] / af     # all 3 arms reach
        print(f"[rc={rc:>4.1f}Å] P(>=2 arms | Anchor) = {p_ge2arms:.3f}  |  P(3 arms | Anchor) = {p_ge3arms:.3f}")

    # Dwell time estimate
    print("\n" + "="*80)
    print("=== Dwell Time (consecutive frames) for ALL-3-arms reach within Anchor ===")
    print("="*80)
    dt = getattr(u.trajectory, "dt", None)
    for rc in REACH_CUTOFFS:
        series = ge2_series[rc] & anchor_series
        # run-length encoding
        lengths = []
        run = 0
        for v in series:
            if v:
                run += 1
            else:
                if run > 0:
                    lengths.append(run)
                    run = 0
        if run > 0:
            lengths.append(run)

        if len(lengths) == 0:
            print(f"[rc={rc:>4.1f}Å] no events")
            continue

        lengths = np.array(lengths, dtype=int)
        if dt is None:
            print(f"[rc={rc:>4.1f}Å] events={len(lengths)}  mean={lengths.mean():.1f} frames  p95={np.percentile(lengths,95):.1f} frames")
        else:
            print(f"[rc={rc:>4.1f}Å] events={len(lengths)}  mean={lengths.mean()*dt:.3f} ps  p95={np.percentile(lengths,95)*dt:.3f} ps")

    # Save arrays for plotting later
    np.savez(
        OUT_NPZ,
        arm_to_target=arm_to_target,
        anchor_arm=anchor_arm,
        anchor_dist=anchor_dist,
        anchor_series=anchor_series.astype(np.int8),
        reach_cutoffs=np.array(REACH_CUTOFFS, dtype=float),
    )
    print(f"\n" + "="*80)
    print(f"Saved: {OUT_NPZ}")
    print("Next: plot histograms / time series from npz if needed.")
    print("="*80 + "\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python analyze_tripod_anchor_reach.py top.psf traj.dcd")
        print("\nExample:")
        print("  python analyze_tripod_anchor_reach.py step5_production.psf production.dcd")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
