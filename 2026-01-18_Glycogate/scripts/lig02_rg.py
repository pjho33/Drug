#!/usr/bin/env python3
"""
(2) LIG Radius of Gyration (Rg) time series
- 뭉침/펼침 정도를 정량화
- 출력: CSV (time_ns, Rg_A)
"""
import os
import numpy as np
from MDAnalysis import Universe

PSF = "/home/pjho3tr/projects/Drug/2026-01-18_Glycogate/data/solution builder/openmm/step3_input.psf"
DCD = "/home/pjho3tr/projects/Drug/2026-01-18_Glycogate/results/md_200ns_rep1/md_rep1_last200ns.dcd"
OUTDIR = "/home/pjho3tr/projects/Drug/2026-01-18_Glycogate/results/md_200ns_rep1"
OUTCSV = os.path.join(OUTDIR, "lig02_rg_last200ns.csv")
LIG_SEL = "resname LIG"

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    u = Universe(PSF, DCD)
    lig = u.select_atoms(LIG_SEL)
    if lig.n_atoms == 0:
        raise SystemExit(f"ERROR: empty ligand selection: {LIG_SEL}")

    rows = []
    for ts in u.trajectory:
        rows.append((ts.time/1000.0, float(lig.radius_of_gyration())))

    np.savetxt(OUTCSV, rows, delimiter=",", header="time_ns,Rg_A", comments="")
    vals = [r[1] for r in rows]
    print("OK: wrote", OUTCSV)
    print(f"frames={len(rows)}  Rg_A: min={min(vals):.3f}  mean={np.mean(vals):.3f}  max={max(vals):.3f}")

if __name__ == "__main__":
    main()
