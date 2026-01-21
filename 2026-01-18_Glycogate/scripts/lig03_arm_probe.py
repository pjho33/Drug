#!/usr/bin/env python3
"""
(3) '팔 3개가 얼마나 뻗었나'를 위한 자동 프로빙
- 정확한 arm별 endpoint 정의에는 원자/서브구조 정보가 필요하지만,
  여기서는 "자동 후보 제시 + 대략적 3-arm span 추정"까지 해준다.

출력:
- lig03_atomnames.txt : LIG 고유 atom name 목록
- lig03_arm_candidates.csv : (time_ns, span1_A, span2_A, span3_A, core_atom_index, tip1_index, tip2_index, tip3_index)
  * core = 프레임0에서 LIG COM에 가장 가까운 atom
  * tips = 프레임0에서 core로부터 가장 먼 원자 3개를 "3개 팔 끝" 후보로 사용
  * 이건 '후보'이며, 실제 arm endpoint로 확정하려면 원자 이름 기반 지정이 최종정답
"""
import os
import numpy as np
from MDAnalysis import Universe

PSF = "/home/pjho3tr/projects/Drug/2026-01-18_Glycogate/data/solution builder/openmm/step3_input.psf"
DCD = "/home/pjho3tr/projects/Drug/2026-01-18_Glycogate/results/md_200ns_rep1/md_rep1_last200ns.dcd"
OUTDIR = "/home/pjho3tr/projects/Drug/2026-01-18_Glycogate/results/md_200ns_rep1"
OUTTXT = os.path.join(OUTDIR, "lig03_atomnames.txt")
OUTCSV = os.path.join(OUTDIR, "lig03_arm_candidates_last200ns.csv")
LIG_SEL = "resname LIG"

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    u = Universe(PSF, DCD)
    lig = u.select_atoms(LIG_SEL)
    if lig.n_atoms == 0:
        raise SystemExit(f"ERROR: empty ligand selection: {LIG_SEL}")

    # (A) atom name 목록 저장
    uniq_names = sorted(set(lig.names))
    with open(OUTTXT, "w") as f:
        f.write(f"LIG atoms: {lig.n_atoms}\n")
        f.write("Unique atom names:\n")
        for n in uniq_names:
            f.write(n + "\n")
    print("OK: wrote", OUTTXT, f"(unique atom names={len(uniq_names)})")

    # (B) 프레임0에서 core/tips 후보 선택
    u.trajectory[0]
    com0 = lig.center_of_mass()
    d0 = np.linalg.norm(lig.positions - com0, axis=1)
    core_i = int(np.argmin(d0))  # COM에 가장 가까운 원자
    core_pos0 = lig.positions[core_i]

    # core로부터 먼 순서 top 3을 tip 후보로
    dc = np.linalg.norm(lig.positions - core_pos0, axis=1)
    tip_indices = np.argsort(dc)[-3:][::-1].astype(int).tolist()
    tip1, tip2, tip3 = tip_indices

    # 후보 원자 정보 출력
    def atom_info(i):
        a = lig.atoms[i]
        return f"idx={i}  name={a.name}  type={a.type}  resid={a.resid}  resname={a.resname}"

    print("Core candidate (frame0):", atom_info(core_i))
    print("Tip candidates (frame0):")
    print("  tip1:", atom_info(tip1))
    print("  tip2:", atom_info(tip2))
    print("  tip3:", atom_info(tip3))

    # (C) time series: core->tip distances
    rows = []
    for ts in u.trajectory:
        t_ns = ts.time / 1000.0
        core_pos = lig.positions[core_i]
        p1 = lig.positions[tip1]
        p2 = lig.positions[tip2]
        p3 = lig.positions[tip3]
        s1 = float(np.linalg.norm(p1 - core_pos))
        s2 = float(np.linalg.norm(p2 - core_pos))
        s3 = float(np.linalg.norm(p3 - core_pos))
        rows.append((t_ns, s1, s2, s3, core_i, tip1, tip2, tip3))

    np.savetxt(OUTCSV, rows, delimiter=",",
               header="time_ns,span1_A,span2_A,span3_A,core_i,tip1_i,tip2_i,tip3_i", comments="")
    spans = np.array([[r[1], r[2], r[3]] for r in rows], dtype=float)
    print("OK: wrote", OUTCSV)
    print(f"frames={len(rows)}  span_A (core->tips):")
    print("  span1 min/mean/max:", float(spans[:,0].min()), float(spans[:,0].mean()), float(spans[:,0].max()))
    print("  span2 min/mean/max:", float(spans[:,1].min()), float(spans[:,1].mean()), float(spans[:,1].max()))
    print("  span3 min/mean/max:", float(spans[:,2].min()), float(spans[:,2].mean()), float(spans[:,2].max()))
    print("NOTE: tips/core are heuristics. If you want true 3-arm endpoints, we will pin atom names from lig03_atomnames.txt.")

if __name__ == "__main__":
    main()
