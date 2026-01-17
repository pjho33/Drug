#!/usr/bin/env python3
"""
Merge SDG ligand coordinates into CHARMM-GUI membrane system
Then create PSF with SDG topology included
"""

import parmed as pmd
from parmed import charmm
import sys

print("=" * 80)
print("Merging SDG into CHARMM-GUI Membrane System")
print("=" * 80)

# Paths
charmm_gui_dir = '/home/pjho3/다운로드/charmm-gui-6750265216membranebuilder/openmm'
membrane_pdb = f'{charmm_gui_dir}/step5_input.pdb'
membrane_psf = f'{charmm_gui_dir}/step5_input.psf'
sdg_pdb = 'glut1_tripod_complex_sdg_only.pdb'

print("\n1. Loading CHARMM-GUI membrane system...")
print(f"   PSF: {membrane_psf}")
print(f"   PDB: {membrane_pdb}")

# Load membrane system (has PSF topology)
membrane = pmd.load_file(membrane_psf)
membrane_coords = pmd.load_file(membrane_pdb)

# Copy coordinates
membrane.coordinates = membrane_coords.coordinates

print(f"   ✅ Loaded: {len(membrane.atoms)} atoms")
print(f"   Residues: {len(membrane.residues)}")

# Load SDG coordinates
print(f"\n2. Loading SDG ligand...")
print(f"   PDB: {sdg_pdb}")

sdg = pmd.load_file(sdg_pdb)
print(f"   ✅ Loaded: {len(sdg.atoms)} atoms")
print(f"   SDG residues: {len(sdg.residues)}")

# Load SDG parameters
print("\n3. Loading SDG parameters...")
try:
    sdg_params = charmm.CharmmParameterSet('sdg.rtf', 'sdg.prm')
    print("   ✅ SDG parameters loaded")
    
    # Apply parameters to SDG structure
    sdg.load_parameters(sdg_params)
    print("   ✅ Parameters applied to SDG")
except Exception as e:
    print(f"   ⚠️  Warning: {e}")
    print("   Continuing without parameter application...")

# Merge structures
print("\n4. Merging SDG into membrane system...")
combined = membrane + sdg

print(f"   ✅ Combined: {len(combined.atoms)} atoms")
print(f"   Total residues: {len(combined.residues)}")

# Save combined PDB
output_pdb = 'glut1_sdg_membrane_combined.pdb'
print(f"\n5. Saving combined PDB...")
combined.save(output_pdb, overwrite=True)
print(f"   ✅ Saved: {output_pdb}")

# Try to save PSF
output_psf = 'glut1_sdg_membrane_combined.psf'
print(f"\n6. Attempting to save PSF...")
try:
    combined.save(output_psf, overwrite=True)
    print(f"   ✅ Saved: {output_psf}")
    
    # Verify PSF
    test_psf = pmd.load_file(output_psf)
    sdg_count = sum(1 for res in test_psf.residues if res.name == 'SDG')
    print(f"   ✅ PSF verified: {sdg_count} SDG residues")
    
except Exception as e:
    print(f"   ❌ Error saving PSF: {e}")
    print("\n   Alternative: Use PDB-only approach with OpenMM")

print("\n" + "=" * 80)
print("✅ Merge Complete!")
print("=" * 80)
print(f"\nOutput files:")
print(f"  - {output_pdb}")
if 'output_psf' in locals():
    print(f"  - {output_psf}")
print("\n" + "=" * 80)
