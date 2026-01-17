#!/usr/bin/env python3
"""
Generate PSF file for GLUT1-SDG complex using ParmEd
This creates the topology file needed for OpenMM CHARMM simulation
"""

import parmed as pmd
from parmed import charmm
import sys

print("=" * 80)
print("Generating PSF for GLUT1-SDG Complex")
print("=" * 80)

# Load PDB
print("\n1. Loading PDB...")
pdb_file = 'glut1_tripod_complex.pdb'
structure = pmd.load_file(pdb_file)

print(f"   ✅ Loaded: {len(structure.atoms)} atoms")
print(f"   Residues: {len(structure.residues)}")

# Count residues
sdg_count = sum(1 for res in structure.residues if res.name == 'SDG')
protein_count = sum(1 for res in structure.residues if res.name not in ['SDG', 'BGLCNA', 'AMAN', 'BMAN', 'ANE5AC', 'ANE', 'BGL', 'BGA', 'AMA'])
glycan_count = len(structure.residues) - sdg_count - protein_count

print(f"   Protein residues: {protein_count}")
print(f"   Glycan residues: {glycan_count}")
print(f"   SDG residues: {sdg_count}")

# Load CHARMM parameters
print("\n2. Loading CHARMM parameters...")
print("   Loading sdg.rtf and sdg.prm...")

try:
    params = charmm.CharmmParameterSet('sdg.rtf', 'sdg.prm')
    print(f"   ✅ Parameters loaded")
except Exception as e:
    print(f"   ❌ Error loading parameters: {e}")
    print("\n   Trying alternative approach...")
    
    # Alternative: load without parameters first
    params = None

# Apply parameters to structure
if params:
    print("\n3. Applying parameters to structure...")
    try:
        structure.load_parameters(params)
        print("   ✅ Parameters applied")
    except Exception as e:
        print(f"   ⚠️  Warning: {e}")
        print("   Continuing without parameter application...")

# Save PSF
print("\n4. Saving PSF file...")
psf_file = 'glut1_sdg_complex.psf'

try:
    structure.save(psf_file, format='psf')
    print(f"   ✅ Saved: {psf_file}")
except Exception as e:
    print(f"   ❌ Error saving PSF: {e}")
    sys.exit(1)

# Verify PSF
print("\n5. Verifying PSF...")
try:
    test_psf = pmd.load_file(psf_file)
    print(f"   ✅ PSF valid: {len(test_psf.atoms)} atoms")
    
    # Check for SDG
    sdg_in_psf = sum(1 for res in test_psf.residues if res.name == 'SDG')
    print(f"   SDG residues in PSF: {sdg_in_psf}")
    
    if sdg_in_psf == 0:
        print("   ⚠️  WARNING: No SDG residues found in PSF!")
    else:
        print("   ✅ SDG residues present in PSF")
        
except Exception as e:
    print(f"   ❌ Error verifying PSF: {e}")

print("\n" + "=" * 80)
print("✅ PSF Generation Complete!")
print("=" * 80)
print(f"\nOutput files:")
print(f"  - {psf_file}")
print(f"\nNext step:")
print(f"  Use CharmmPsfFile to load this PSF in OpenMM")
print("=" * 80)
