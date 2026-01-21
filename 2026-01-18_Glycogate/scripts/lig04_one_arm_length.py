#!/usr/bin/env python3
"""
1-arm 리간드용: core(heavy atom) -> tip(heavy atom) 최대거리 time series

- core: frame0에서 LIG heavy-atom COM에 가장 가까운 원자
- tip : frame0에서 core로부터 가장 먼 heavy atom 1개 (팔 끝 후보)
- 출력: CSV (time_ns, arm_A)
- 4nm/5nm(40/50Å) 이상 비율 출력
"""
import os
import numpy as np
from MDAnalysis import Universe
from MDAnalysis.lib.distances import calc_bonds

PSF = "/home/pjho3tr/projects/Drug/2026-01-18_Glycogate/data/solution builder/openmm/step3_input.psf"
DCD = "/home/pjho3tr/projects/Drug/2026-01-18_Glycogate/results/md_200ns_rep1/md_rep1_last200ns.dcd"
OUTDIR = "/home/pjho3tr/projects/Drug/2026-01-18_Glycogate/results/md_200ns_rep1"
OUTCSV = os.path.join(OUTDIR, "lig04_one_arm_length_last200ns.csv")

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

    # choose core/tip on frame0
    u.trajectory[0]
    com = lig.center_of_mass()
    d_com = np.linalg.norm(lig.positions - com, axis=1)
    core_local = int(np.argmin(d_com))
    core = lig.atoms[core_local]

    d_core = np.linalg.norm(lig.positions - core.position, axis=1)
    tip_local = int(np.argmax(d_core))
    tip = lig.atoms[tip_local]

    print("Selected CORE:", info(core))
    print("Selected TIP :", info(tip))

    core_ag = u.atoms[core.ix:core.ix+1]
    tip_ag  = u.atoms[tip.ix:tip.ix+1]

    rows=[]
    vals=[]
    for ts in u.trajectory:
        dist = float(calc_bonds(core_ag.positions, tip_ag.positions, box=ts.dimensions)[0])
        t_ns = ts.time/1000.0
        rows.append((t_ns, dist))
        vals.append(dist)

    with open(OUTCSV, "w") as f:
        f.write("time_ns,arm_A\n")
        for t,d in rows:
            f.write(f"{t:.6f},{d:.6f}\n")

    vals=np.array(vals)
    print("OK: wrote", OUTCSV)
    print(f"frames={len(vals)} arm_A min/mean/max = {vals.min():.3f} / {vals.mean():.3f} / {vals.max():.3f}")
    print(f"P(arm >= 40Å / 4nm) [%] = {(vals>=THR40).mean()*100:.2f}")
    print(f"P(arm >= 50Å / 5nm) [%] = {(vals>=THR50).mean()*100:.2f}")

if __name__ == "__main__":
    main()
