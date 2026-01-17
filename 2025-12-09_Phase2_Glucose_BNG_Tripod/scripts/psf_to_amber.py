#!/usr/bin/env python3
"""
Convert PSF + toppar to Amber prmtop/rst7
Using ParmEd to convert CHARMM topology to Amber format
"""

import argparse
import glob
import parmed as pmd
from parmed.charmm import CharmmPsfFile, CharmmParameterSet
import os

def collect_toppar_files(toppar_path: str, ligand_dir: str = None) -> list:
    # If it's a single file, return it
    if os.path.isfile(toppar_path):
        return [toppar_path]
    
    # If it's a directory, search for toppar files
    files = []
    if os.path.isdir(toppar_path):
        pats = ["*.str", "*.prm", "*.rtf", "*.par"]
        for pat in pats:
            files.extend(sorted(glob.glob(os.path.join(toppar_path, pat))))
    
    # Add ligand-specific toppar files if provided
    if ligand_dir and os.path.isdir(ligand_dir):
        pats = ["*.str", "*.prm", "*.rtf", "*.par"]
        for pat in pats:
            files.extend(sorted(glob.glob(os.path.join(ligand_dir, pat))))
    
    return files

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--psf", required=True, help="Dry PSF file")
    ap.add_argument("--pdb", required=True, help="Dry PDB (protein+glycan+ligand)")
    ap.add_argument("--toppar", required=True, help="CHARMM-GUI toppar file or directory")
    ap.add_argument("--ligand-toppar", help="Ligand-specific toppar directory (e.g., sdg/)")
    ap.add_argument("--out", default="complex", help="Output prefix (default: complex)")
    args = ap.parse_args()

    print("="*70)
    print("PSF → Amber prmtop Conversion")
    print("="*70)

    psf_file = args.psf
    pdb_file = args.pdb
    toppar_path = args.toppar
    ligand_toppar = args.ligand_toppar

    print(f"\n📂 Input files:")
    print(f"  PSF: {psf_file}")
    print(f"  PDB: {pdb_file}")
    print(f"  Toppar: {toppar_path}")
    if ligand_toppar:
        print(f"  Ligand toppar: {ligand_toppar}")

    # Find toppar files
    print(f"\n🔍 Loading toppar files...")
    toppar_files = collect_toppar_files(toppar_path, ligand_toppar)

    if not toppar_files:
        raise SystemExit(f"❌ No toppar files found at {toppar_path} (need .str/.prm/.rtf etc.)")

    print(f"✅ Found {len(toppar_files)} toppar file(s)")

    # Load CHARMM parameters
    print(f"\n📖 Loading CHARMM parameters...")
    params = CharmmParameterSet(*toppar_files)
    print(f"✅ Parameters loaded")

    # Load PSF
    print(f"\n📖 Loading PSF...")
    psf = CharmmPsfFile(psf_file)
    psf.load_parameters(params)
    print(f"✅ PSF loaded: {len(psf.atoms):,} atoms")

    # Load coordinates
    print(f"\n📖 Loading coordinates...")
    coords = pmd.load_file(pdb_file)
    
    # Check atom count match
    if len(psf.atoms) != len(coords.atoms):
        raise SystemExit(
            f"❌ Atom count mismatch:\n"
            f"  PSF natom={len(psf.atoms)}\n"
            f"  PDB natom={len(coords.atoms)}\n"
            f"Dry PDB must correspond to the SAME atom ordering as PSF."
        )
    
    psf.coordinates = coords.coordinates
    print(f"✅ Coordinates loaded: {len(coords.atoms):,} atoms")

    # Save as Amber format
    out_parm7 = f"{args.out}.parm7"
    out_rst7 = f"{args.out}.rst7"
    
    print(f"\n💾 Saving Amber files...")
    psf.save(out_parm7, overwrite=True)
    psf.save(out_rst7, overwrite=True)
    
    print(f"✅ Wrote {out_parm7} / {out_rst7}")
    print(f"\n{'='*70}")

if __name__ == "__main__":
    main()
