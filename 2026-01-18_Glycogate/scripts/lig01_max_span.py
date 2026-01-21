#!/usr/bin/env python3
"""
(1) LIG 최대 internal span (프레임별 LIG 내부 원자쌍 최대거리)
- "팔을 얼마나 뻗치냐"의 가장 직관적인 proxy
- 출력: CSV (time_ns, max_internal_distance_A)
"""
import os
import numpy as np
from MDAnalysis import Universe
from MDAnalysis.lib.distances import distance_array

PSF = "/home/pjho3tr/projects/Drug/2026-01-18_Glycogate/data/solution builder/openmm/step3_input.psf"
DCD = "/home/pjho3tr/projects/Drug/2026-01-18_Glycogate/results/md_200ns_rep1/md_rep1_last200ns.dcd"
OUTDIR = "/home/pjho3tr/projects/Drug/2026-01-18_Glycogate/results/md_200ns_rep1"
OUTCSV = os.path.join(OUTDIR, "lig01_max_span_last200ns.csv")
LIG_SEL = "resname LIG"

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    u = Universe(PSF, DCD)
    lig = u.select_atoms(LIG_SEL)
    if lig.n_atoms == 0:
        raise SystemExit(f"ERROR: empty ligand selection: {LIG_SEL}")

    rows = []
    for ts in u.trajectory:
        # N x N distance matrix, take max
        D = distance_array(lig.positions, lig.positions, box=ts.dimensions)
        maxd = float(D.max())
        rows.append((ts.time/1000.0, maxd))

    np.savetxt(OUTCSV, rows, delimiter=",",
               header="time_ns,max_internal_distance_A", comments="")
    vals = [r[1] for r in rows]
    print("OK: wrote", OUTCSV)
    print(f"frames={len(rows)}  max_span_A: min={min(vals):.3f}  mean={np.mean(vals):.3f}  max={max(vals):.3f}")

if __name__ == "__main__":
    main()
