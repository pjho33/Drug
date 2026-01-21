#!/usr/bin/env python3
"""
lig04_true_arm_lengths.py

목적:
- LIG(=TRIS-PEG24-L-glucose, 3-arm)의 '코어 + 말단 3개'를 heavy atom 기준으로 자동 선택
- 시간에 따른 arm length (core -> tip1/tip2/tip3) 계산
- 4-5 nm(=40-50 Å) 기준 초과 비율 및 동시 초과 비율 출력

입력(고정 경로, 필요시 직접 수정):
- PSF, DCD

출력:
- lig04_true_arm_lengths_last200ns.csv
- 콘솔 요약(선택된 core/tips atom info + 통계 + threshold 통과율)
"""

import os
import numpy as np
from MDAnalysis import Universe
from MDAnalysis.lib.distances import calc_bonds

# ====== 사용자 환경 경로(필요시 수정) ======
PSF = "/home/pjho3tr/projects/Drug/2026-01-18_Glycogate/data/solution builder/openmm/step3_input.psf"
DCD = "/home/pjho3tr/projects/Drug/2026-01-18_Glycogate/results/md_200ns_rep1/md_rep1_last200ns.dcd"
OUTDIR = "/home/pjho3tr/projects/Drug/2026-01-18_Glycogate/results/md_200ns_rep1"
OUTCSV = os.path.join(OUTDIR, "lig04_true_arm_lengths_last200ns.csv")

# LIG selection (heavy atoms only)
LIG_HEAVY_SEL = "resname LIG and not name H*"

# tip selection heuristic parameters
TOPN_CANDIDATES = 200     # core에서 먼 원자 상위 N개 중에서 tip 3개 선택
MIN_TIP_SEP_A = 15.0      # tip들끼리 서로 최소 이 정도는 떨어져야 "다른 arm" 가능성↑
MIN_CORE_DIST_A = 10.0    # tip 후보는 core로부터 최소 이 정도는 떨어져야 함(너무 가까운 것 배제)

# threshold checks (nm -> Å)
THR1_A = 40.0  # 4 nm
THR2_A = 50.0  # 5 nm


def atom_info(atom):
    return f"ix={atom.ix} name={atom.name} type={atom.type} resid={atom.resid} resname={atom.resname}"


def pick_core_and_tips(lig_heavy):
    """
    frame0 기준:
    - core: COM에 가장 가까운 heavy atom
    - tips: core로부터 멀고, tip끼리 떨어진(>= MIN_TIP_SEP_A) heavy atom 3개
    """
    # core = closest-to-COM heavy atom
    com = lig_heavy.center_of_mass()
    d_com = np.linalg.norm(lig_heavy.positions - com, axis=1)
    core_local_i = int(np.argmin(d_com))
    core_atom = lig_heavy.atoms[core_local_i]

    core_pos = core_atom.position
    d_core = np.linalg.norm(lig_heavy.positions - core_pos, axis=1)

    # 후보: core로부터 충분히 먼 heavy atom들 중, 먼 순으로 정렬
    cand_local = np.where(d_core >= MIN_CORE_DIST_A)[0]
    if cand_local.size == 0:
        raise RuntimeError("No candidates far from core. Try lowering MIN_CORE_DIST_A.")

    # 거리 큰 순
    cand_local = cand_local[np.argsort(d_core[cand_local])[::-1]]
    cand_local = cand_local[:min(TOPN_CANDIDATES, len(cand_local))]

    # greedy tip selection with separation constraint
    tips_local = []
    for idx in cand_local:
        if not tips_local:
            tips_local.append(int(idx))
            if len(tips_local) == 3:
                break
            continue

        pos = lig_heavy.positions[idx]
        ok = True
        for t in tips_local:
            if np.linalg.norm(pos - lig_heavy.positions[t]) < MIN_TIP_SEP_A:
                ok = False
                break
        if ok:
            tips_local.append(int(idx))
            if len(tips_local) == 3:
                break

    # 만약 3개 못 고르면 separation 완화하며 재시도
    if len(tips_local) < 3:
        for relaxed in [12.0, 10.0, 8.0, 6.0]:
            tips_local = []
            for idx in cand_local:
                if not tips_local:
                    tips_local.append(int(idx))
                    if len(tips_local) == 3:
                        break
                    continue
                pos = lig_heavy.positions[idx]
                ok = True
                for t in tips_local:
                    if np.linalg.norm(pos - lig_heavy.positions[t]) < relaxed:
                        ok = False
                        break
                if ok:
                    tips_local.append(int(idx))
                    if len(tips_local) == 3:
                        break
            if len(tips_local) == 3:
                break

    if len(tips_local) != 3:
        raise RuntimeError(
            "Failed to pick 3 separated tips. "
            "Increase TOPN_CANDIDATES or lower MIN_TIP_SEP_A / MIN_CORE_DIST_A."
        )

    tip_atoms = [lig_heavy.atoms[i] for i in tips_local]
    return core_atom, tip_atoms


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    u = Universe(PSF, DCD)
    lig_heavy = u.select_atoms(LIG_HEAVY_SEL)
    if lig_heavy.n_atoms == 0:
        raise SystemExit(f"ERROR: empty ligand selection: {LIG_HEAVY_SEL}")

    # frame0에서 core/tips 선택
    u.trajectory[0]
    core_atom, tip_atoms = pick_core_and_tips(lig_heavy)

    print("Selected CORE (heavy atom, frame0):", atom_info(core_atom))
    for i, t in enumerate(tip_atoms, 1):
        print(f"Selected TIP{i}  (heavy atom, frame0):", atom_info(t))

    # core/tip local indices within lig_heavy (for fast indexing)
    # Use atom.ix (global index) then map to local positions by selection
    # Easiest: just use atom indices directly and calc_bonds with arrays of positions each frame
    core_ix = core_atom.ix
    tip_ixs = [t.ix for t in tip_atoms]

    # Build AtomGroups for stable referencing
    core_ag = u.atoms[core_ix:core_ix+1]
    tip_ags = [u.atoms[i:i+1] for i in tip_ixs]

    rows = []
    arm_lengths = []  # list of (a1,a2,a3)

    for ts in u.trajectory:
        box = ts.dimensions
        cpos = core_ag.positions
        a = []
        for tag in tip_ags:
            # calc_bonds: returns array of distances
            dist = float(calc_bonds(cpos, tag.positions, box=box)[0])
            a.append(dist)
        t_ns = ts.time / 1000.0
        rows.append((t_ns, a[0], a[1], a[2]))
        arm_lengths.append(a)

    arm_lengths = np.array(arm_lengths, dtype=float)

    # save CSV
    with open(OUTCSV, "w") as f:
        f.write("time_ns,arm1_A,arm2_A,arm3_A\n")
        for r in rows:
            f.write(f"{r[0]:.6f},{r[1]:.6f},{r[2]:.6f},{r[3]:.6f}\n")

    print("OK: wrote", OUTCSV)
    print("frames:", len(rows))
    for k in range(3):
        v = arm_lengths[:, k]
        print(f"arm{k+1} Å  min/mean/max = {v.min():.3f} / {v.mean():.3f} / {v.max():.3f}")

    # threshold stats
    thr1 = (arm_lengths >= THR1_A)
    thr2 = (arm_lengths >= THR2_A)

    p_arm_ge_4nm = thr1.mean(axis=0) * 100.0
    p_arm_ge_5nm = thr2.mean(axis=0) * 100.0

    print(f"\nThreshold pass rates (per arm):")
    print(f"  P(arm >= 4.0 nm / 40Å) [%] = {p_arm_ge_4nm[0]:.1f}, {p_arm_ge_4nm[1]:.1f}, {p_arm_ge_4nm[2]:.1f}")
    print(f"  P(arm >= 5.0 nm / 50Å) [%] = {p_arm_ge_5nm[0]:.1f}, {p_arm_ge_5nm[1]:.1f}, {p_arm_ge_5nm[2]:.1f}")

    # simultaneous: at least 2 arms exceed threshold
    p_atleast2_4nm = (thr1.sum(axis=1) >= 2).mean() * 100.0
    p_atleast2_5nm = (thr2.sum(axis=1) >= 2).mean() * 100.0
    p_all3_4nm = (thr1.sum(axis=1) == 3).mean() * 100.0
    p_all3_5nm = (thr2.sum(axis=1) == 3).mean() * 100.0

    print(f"\nSimultaneous reach (same frame):")
    print(f"  P(at least 2 arms >= 4.0 nm) [%] = {p_atleast2_4nm:.1f}")
    print(f"  P(all 3 arms   >= 4.0 nm) [%] = {p_all3_4nm:.1f}")
    print(f"  P(at least 2 arms >= 5.0 nm) [%] = {p_atleast2_5nm:.1f}")
    print(f"  P(all 3 arms   >= 5.0 nm) [%] = {p_all3_5nm:.1f}")

    print("\nNOTE:")
    print("- core/tips are chosen by geometric heuristics on frame0 (heavy atoms only).")
    print("- If you want chemically exact endpoints (TRIS core atom + each terminal glucose atom),")
    print("  we can lock them by atom names/indices once we identify them from lig03_atomnames.txt.")

if __name__ == "__main__":
    main()
