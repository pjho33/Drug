#!/usr/bin/env python3
"""
선택자 확인용 스크립트: Glucose/리간드 resid 찾기
Tripod의 3개 말단 Glucose와 Target residue 범위를 확인
"""

import sys
import MDAnalysis as mda
from collections import Counter

def main(top, traj=None):
    u = mda.Universe(top, traj) if traj else mda.Universe(top)

    # 전체 residue 통계
    resnames = [r.resname for r in u.residues]
    c = Counter(resnames)
    print("\n" + "="*80)
    print("[Top 30 resnames by count]")
    print("="*80)
    for k, v in c.most_common(30):
        print(f"{k:>8s} : {v:>6d}")

    # Glucose 후보를 넓게 탐색 (이름이 다를 수 있어서)
    candidates = []
    for r in u.residues:
        rn = r.resname.upper()
        if any(key in rn for key in ["GLC", "BGLC", "AGLC", "GLU", "BGC", "SUG", "GL"]):
            candidates.append((r.resid, r.resname, r.n_atoms))

    print("\n" + "="*80)
    print("[Residues that look like sugars (heuristic)]")
    print("="*80)
    if not candidates:
        print("  (none found by heuristic) -> you must search by ligand/segment name manually")
    else:
        for resid, resname, natoms in candidates[:200]:
            print(f"  resid {resid:<6d} resname {resname:<8s} natoms {natoms}")

    # 리간드처럼 원자 수가 작은 residue 후보(5~400 atoms)를 나열
    small = [(r.resid, r.resname, r.n_atoms) for r in u.residues if 5 <= r.n_atoms <= 400]
    print("\n" + "="*80)
    print("[Small residues (5~400 atoms) - first 50]")
    print("="*80)
    for resid, resname, natoms in small[:50]:
        print(f"  resid {resid:<6d} resname {resname:<8s} natoms {natoms}")

    # Protein residues 범위 확인
    protein = u.select_atoms("protein")
    if len(protein) > 0:
        prot_resids = sorted(set([r.resid for r in protein.residues]))
        print("\n" + "="*80)
        print("[Protein residue range]")
        print("="*80)
        print(f"  Total protein residues: {len(prot_resids)}")
        print(f"  Range: {min(prot_resids)} - {max(prot_resids)}")
        print(f"  First 20: {prot_resids[:20]}")
        print(f"  Last 20: {prot_resids[-20:]}")

    print("\n" + "="*80)
    print("TIP: If you know ligand resname (e.g., SDG, LIG), run selection like:")
    print("  u.select_atoms('resname SDG') and inspect residues/resids.")
    print("\nFor Tripod analysis, you need:")
    print("  1. Three Glucose resids (말단 3개)")
    print("  2. Target residue range (vestibule 주변)")
    print("="*80 + "\n")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python 00_inspect_selections.py top.psf [traj.dcd]")
        print("\nExample:")
        print("  python 00_inspect_selections.py step5_production.psf")
        print("  python 00_inspect_selections.py step5_production.psf production.dcd")
        sys.exit(1)
    top = sys.argv[1]
    traj = sys.argv[2] if len(sys.argv) >= 3 else None
    main(top, traj)
