#!/usr/bin/env python3
"""
Create dry PSF (protein + glycan + ligand only)
Remove water, ions, lipids from PSF while maintaining atom ordering
"""

import argparse
import parmed as pmd
import os

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--psf", required=True, help="Original PSF file")
    ap.add_argument("--pdb", required=True, help="Full system PDB matching PSF atom order")
    ap.add_argument("--out", default="dry_complex", help="Output prefix")
    ap.add_argument("--lig", default="SDG", help="Ligand residue name")
    ap.add_argument("--glycan", nargs="*", default=["NAG","BMA","MAN","GAL","FUC","SIA","NEU5AC","NDG","BGLC","A2G","A2M"], 
                    help="Glycan residue names")
    args = ap.parse_args()

    print("="*70)
    print("Creating Dry PSF (protein + glycan + ligand)")
    print("="*70)

    # Load PSF + coordinates
    print(f"\n📂 Loading PSF and PDB...")
    st = pmd.load_file(args.psf)
    pdb = pmd.load_file(args.pdb)
    st.coordinates = pdb.coordinates
    print(f"✅ Loaded: {len(st.atoms):,} atoms")

    glycan_set = set(args.glycan)
    lig = args.lig

    # Exclude list: water, ions, lipids
    exclude_resnames = set([
        "TIP3", "TIP3P", "WAT", "HOH",
        "NA", "CL", "K", "CA", "MG", "POT", "CLA", "SOD",
        "POPC", "POPE", "POPS", "POPG", "CHL1", "CHOL",
        "DLPC", "DPPC", "DOPC", "DOPE", "DOPS", "DPPE", "DPPG"
    ])

    # Build mask: keep protein, glycan, ligand; exclude water/ions/lipids
    print(f"\n🔍 Building atom mask...")
    mask = []
    stats = {"protein": 0, "glycan": 0, "ligand": 0, "excluded": 0}
    
    for a in st.atoms:
        rn = a.residue.name.strip()
        
        # Ligand
        if rn == lig:
            mask.append(True)
            stats["ligand"] += 1
            continue
        
        # Glycan
        if rn in glycan_set:
            mask.append(True)
            stats["glycan"] += 1
            continue
        
        # Excluded (water/ions/lipids)
        if rn in exclude_resnames:
            mask.append(False)
            stats["excluded"] += 1
            continue
        
        # Everything else is protein
        mask.append(True)
        stats["protein"] += 1

    print(f"\n📊 Atom statistics:")
    print(f"  Protein: {stats['protein']:,}")
    print(f"  Glycan: {stats['glycan']:,}")
    print(f"  Ligand ({lig}): {stats['ligand']:,}")
    print(f"  Excluded: {stats['excluded']:,}")
    print(f"  Total kept: {sum(mask):,}")

    # Create dry structure
    print(f"\n🔧 Creating dry structure...")
    dry = st[mask]

    # Save dry PSF and PDB
    output_psf = f"{args.out}.psf"
    output_pdb = f"{args.out}.pdb"
    
    print(f"\n💾 Saving dry files...")
    dry.save(output_psf, overwrite=True)
    dry.save(output_pdb, overwrite=True)
    
    print(f"✅ Wrote:")
    print(f"  {output_psf}")
    print(f"  {output_pdb}")
    print(f"  Atoms: {len(dry.atoms):,}")

    print(f"\n{'='*70}")
    print("Dry PSF creation complete!")
    print(f"{'='*70}")

if __name__ == "__main__":
    main()
