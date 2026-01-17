#!/usr/bin/env python3
"""
Rename files only (content already changed by user)
"""

import os
import shutil

print("=" * 80)
print("Renaming files: trp → sdg")
print("=" * 80)

# 1. Backup and rename trp.rtf → sdg.rtf
print("\n1. trp.rtf → sdg.rtf")
if os.path.exists('trp.rtf'):
    shutil.copy('trp.rtf', 'trp_backup.rtf')
    shutil.move('trp.rtf', 'sdg.rtf')
    print("   ✅ Renamed to sdg.rtf")
    print("   ✅ Backup: trp_backup.rtf")
else:
    print("   ⚠️  trp.rtf not found")

# 2. Backup and rename trp.prm → sdg.prm
print("\n2. trp.prm → sdg.prm")
if os.path.exists('trp.prm'):
    shutil.copy('trp.prm', 'trp_backup.prm')
    shutil.move('trp.prm', 'sdg.prm')
    print("   ✅ Renamed to sdg.prm")
    print("   ✅ Backup: trp_backup.prm")
else:
    print("   ⚠️  trp.prm not found")

# 3. Backup PDB (no rename needed)
print("\n3. Backing up glut1_tripod_complex.pdb")
if os.path.exists('glut1_tripod_complex.pdb'):
    shutil.copy('glut1_tripod_complex.pdb', 'glut1_tripod_complex_backup.pdb')
    print("   ✅ Backup: glut1_tripod_complex_backup.pdb")
else:
    print("   ⚠️  glut1_tripod_complex.pdb not found")

print("\n" + "=" * 80)
print("✅ COMPLETE!")
print("=" * 80)
print("\nFiles:")
print("  - sdg.rtf")
print("  - sdg.prm")
print("  - glut1_tripod_complex.pdb (content already updated by user)")
print("\n" + "=" * 80)
