#!/usr/bin/env python3
"""
Rename TRP to SDG in all files to avoid confusion with Tryptophan
"""

import os
import shutil

print("=" * 80)
print("Renaming TRP → SDG (to avoid confusion with Tryptophan)")
print("=" * 80)

# 1. Rename PDB file and replace TRP with SDG
print("\n1. Processing glut1_tripod_complex.pdb...")
with open('glut1_tripod_complex.pdb', 'r') as f:
    content = f.read()

# Replace TRP residue name with SDG
content_sdg = content.replace(' TRP ', ' SDG ')

# Backup original
shutil.copy('glut1_tripod_complex.pdb', 'glut1_tripod_complex_trp_backup.pdb')
print("   ✅ Backup: glut1_tripod_complex_trp_backup.pdb")

# Write new version
with open('glut1_tripod_complex.pdb', 'w') as f:
    f.write(content_sdg)
print("   ✅ Updated: glut1_tripod_complex.pdb (TRP → SDG)")

# Count changes
trp_count = content.count(' TRP ')
print(f"   Changed {trp_count} residues")

# 2. Process trp.rtf → sdg.rtf
print("\n2. Processing trp.rtf → sdg.rtf...")
with open('trp.rtf', 'r') as f:
    rtf_content = f.read()

# Replace RESI TRP with RESI SDG
rtf_sdg = rtf_content.replace('RESI TRP', 'RESI SDG')

with open('sdg.rtf', 'w') as f:
    f.write(rtf_sdg)
print("   ✅ Created: sdg.rtf")

# Backup and remove old
shutil.copy('trp.rtf', 'trp_backup.rtf')
os.remove('trp.rtf')
print("   ✅ Backup: trp_backup.rtf")
print("   ✅ Removed: trp.rtf")

# 3. Rename trp.prm → sdg.prm
print("\n3. Renaming trp.prm → sdg.prm...")
shutil.copy('trp.prm', 'trp_backup.prm')
shutil.move('trp.prm', 'sdg.prm')
print("   ✅ Backup: trp_backup.prm")
print("   ✅ Renamed: sdg.prm")

print("\n" + "=" * 80)
print("✅ COMPLETE!")
print("=" * 80)
print("\nNew files:")
print("  - glut1_tripod_complex.pdb (with SDG residues)")
print("  - sdg.rtf (topology)")
print("  - sdg.prm (parameters)")
print("\nBackup files:")
print("  - glut1_tripod_complex_trp_backup.pdb")
print("  - trp_backup.rtf")
print("  - trp_backup.prm")
print("\n" + "=" * 80)
