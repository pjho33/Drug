#!/usr/bin/env python3
"""
Add mbondi2 radii to Amber prmtop for GB calculations
"""

import argparse
import parmed as pmd
from parmed.tools.actions import changeRadii

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="Input parm7 file")
    ap.add_argument("--rst", required=True, help="Input rst7 file")
    ap.add_argument("--out", required=True, help="Output prefix")
    args = ap.parse_args()

    print("="*70)
    print("Adding mbondi2 radii to prmtop")
    print("="*70)

    print(f"\n📂 Loading files...")
    parm = pmd.load_file(args.inp, args.rst)
    print(f"✅ Loaded: {len(parm.atoms):,} atoms")

    print(f"\n🔧 Applying mbondi2 radii...")
    changeRadii(parm, "mbondi2").execute()
    print(f"✅ Radii applied")

    print(f"\n💾 Saving files...")
    parm.save(f"{args.out}.parm7", overwrite=True)
    parm.save(f"{args.out}.rst7", overwrite=True)
    
    print(f"✅ Wrote:")
    print(f"  {args.out}.parm7")
    print(f"  {args.out}.rst7")
    print(f"\n{'='*70}")

if __name__ == "__main__":
    main()
