#!/usr/bin/env python3
"""
Convert dry PDB to Amber prmtop using tleap
For systems where PSF is not available or atom ordering differs
"""

import argparse
import os
import subprocess

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True, help="Dry PDB file")
    ap.add_argument("--out", default="complex", help="Output prefix")
    args = ap.parse_args()

    print("="*70)
    print("PDB → Amber prmtop via tleap")
    print("="*70)

    # Create tleap input script
    tleap_script = f"""
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p

# Load PDB
mol = loadpdb {args.pdb}

# Save Amber files
saveamberparm mol {args.out}.parm7 {args.out}.rst7
quit
"""
    
    script_file = "tleap_dry.in"
    with open(script_file, "w") as f:
        f.write(tleap_script)
    
    print(f"\n📝 Created tleap script: {script_file}")
    print(f"\n🔧 Running tleap...")
    
    result = subprocess.run(
        ["tleap", "-f", script_file],
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        print(f"❌ tleap failed:")
        print(result.stderr)
        return 1
    
    print(f"✅ tleap completed")
    print(f"\n💾 Output files:")
    print(f"  {args.out}.parm7")
    print(f"  {args.out}.rst7")
    
    # Verify
    if os.path.exists(f"{args.out}.parm7"):
        import parmed as pmd
        parm = pmd.load_file(f"{args.out}.parm7")
        print(f"\n✅ Verification: {len(parm.atoms):,} atoms")
    
    print(f"\n{'='*70}")

if __name__ == "__main__":
    main()
