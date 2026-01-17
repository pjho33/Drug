#!/usr/bin/env python3
"""
Reorganize Drug_local_20160117 into structured folders
- Phase 2: Glucose, BNG, Tripod (replica 1-3)
- Phase 3: SDG-Glycated GLUT1 MD simulation
"""

import os
import shutil
from pathlib import Path

print("="*70)
print("Reorganizing Drug_local_20160117 into structured folders")
print("="*70)

# Paths
SOURCE_DIR = Path("/home/pjho3tr/projects/Drug_local_20160117")
TARGET_BASE = Path("/home/pjho3tr/projects/Drug")

# Create Phase 2 structure
PHASE2_DIR = TARGET_BASE / "20241209_Phase2_Glucose_BNG_Tripod"
PHASE2_SCRIPTS = PHASE2_DIR / "scripts"
PHASE2_RESULTS = PHASE2_DIR / "results"
PHASE2_ANALYSIS = PHASE2_DIR / "analysis"
PHASE2_STRUCTURES = PHASE2_DIR / "structures"

# Create Phase 3 structure
PHASE3_DIR = TARGET_BASE / "20260106_Phase3_SDG_Glycated_GLUT1"
PHASE3_SCRIPTS = PHASE3_DIR / "scripts"
PHASE3_MD = PHASE3_DIR / "md_simulation"
PHASE3_ANALYSIS = PHASE3_DIR / "analysis"
PHASE3_MMPBSA = PHASE3_DIR / "mmpbsa"

print("\n📁 Creating directory structure...")

# Create Phase 2 directories
for dir_path in [PHASE2_SCRIPTS, PHASE2_RESULTS, PHASE2_ANALYSIS, PHASE2_STRUCTURES]:
    dir_path.mkdir(parents=True, exist_ok=True)
    print(f"  ✅ {dir_path}")

# Create Phase 3 directories
for dir_path in [PHASE3_SCRIPTS, PHASE3_MD, PHASE3_ANALYSIS, PHASE3_MMPBSA]:
    dir_path.mkdir(parents=True, exist_ok=True)
    print(f"  ✅ {dir_path}")

print("\n📦 Phase 2: Copying data...")

# Copy Phase 2 results
if (SOURCE_DIR / "results").exists():
    print("  Copying results folder...")
    for item in (SOURCE_DIR / "results").iterdir():
        dest = PHASE2_RESULTS / item.name
        if item.is_dir():
            if not dest.exists():
                shutil.copytree(item, dest)
                print(f"    ✅ {item.name}")
        else:
            shutil.copy2(item, dest)

# Copy Phase 2 scripts
if (SOURCE_DIR / "scripts").exists():
    print("  Copying scripts folder...")
    for item in (SOURCE_DIR / "scripts").iterdir():
        dest = PHASE2_SCRIPTS / item.name
        if item.is_dir():
            if not dest.exists():
                shutil.copytree(item, dest)
        else:
            shutil.copy2(item, dest)
    print(f"    ✅ scripts")

# Copy Phase 2 analysis
if (SOURCE_DIR / "analysis").exists():
    print("  Copying analysis folder...")
    for item in (SOURCE_DIR / "analysis").iterdir():
        dest = PHASE2_ANALYSIS / item.name
        if item.is_dir():
            if not dest.exists():
                shutil.copytree(item, dest)
        else:
            shutil.copy2(item, dest)
    print(f"    ✅ analysis")

# Copy Phase 2 structures
if (SOURCE_DIR / "structures").exists():
    print("  Copying structures folder...")
    for item in (SOURCE_DIR / "structures").iterdir():
        dest = PHASE2_STRUCTURES / item.name
        if item.is_dir():
            if not dest.exists():
                shutil.copytree(item, dest)
        else:
            shutil.copy2(item, dest)
    print(f"    ✅ structures")

# Copy ligand SMILES files
print("  Copying ligand SMILES files...")
for smi_file in SOURCE_DIR.glob("*.smi"):
    shutil.copy2(smi_file, PHASE2_STRUCTURES / smi_file.name)
    print(f"    ✅ {smi_file.name}")

for sdf_file in SOURCE_DIR.glob("*.sdf"):
    shutil.copy2(sdf_file, PHASE2_STRUCTURES / sdf_file.name)
    print(f"    ✅ {sdf_file.name}")

# Copy pipeline scripts
for script in ["run_pipeline.sh", "run_manual_pipeline.sh", "run_pipeline_now.sh", "run_tris_peg12.sh"]:
    if (SOURCE_DIR / script).exists():
        shutil.copy2(SOURCE_DIR / script, PHASE2_SCRIPTS / script)
        print(f"    ✅ {script}")

print("\n📦 Phase 3: Organizing SDG-Glycated GLUT1 data...")

# Copy Phase 3 archive data
if (SOURCE_DIR / "archive").exists():
    print("  Copying archive data...")
    for item in (SOURCE_DIR / "archive").iterdir():
        if "phase3" in item.name.lower():
            dest = PHASE3_MD / item.name
            if item.is_dir() and not dest.exists():
                shutil.copytree(item, dest)
                print(f"    ✅ {item.name}")

# Copy Phase 3 final data
if (SOURCE_DIR / "phase3_final").exists():
    print("  Copying phase3_final...")
    dest = PHASE3_MD / "phase3_final"
    if not dest.exists():
        shutil.copytree(SOURCE_DIR / "phase3_final", dest)
        print(f"    ✅ phase3_final")

# Copy MMPBSA data
if (SOURCE_DIR / "mmpbsa").exists():
    print("  Copying mmpbsa...")
    for item in (SOURCE_DIR / "mmpbsa").iterdir():
        dest = PHASE3_MMPBSA / item.name
        if item.is_dir():
            if not dest.exists():
                shutil.copytree(item, dest)
        else:
            shutil.copy2(item, dest)
    print(f"    ✅ mmpbsa")

# Copy Phase 3 MD analysis
if (SOURCE_DIR / "phase3_md_analysis").exists():
    print("  Copying phase3_md_analysis...")
    for item in (SOURCE_DIR / "phase3_md_analysis").iterdir():
        dest = PHASE3_ANALYSIS / item.name
        if item.is_dir():
            if not dest.exists():
                shutil.copytree(item, dest)
        else:
            shutil.copy2(item, dest)

# Copy current Drug folder scripts to Phase 3
CURRENT_SCRIPTS = Path("/home/pjho3tr/projects/Drug/scripts")
if CURRENT_SCRIPTS.exists():
    print("  Copying current Drug/scripts to Phase 3...")
    for item in CURRENT_SCRIPTS.iterdir():
        if item.name.startswith(('.', '__')):
            continue
        dest = PHASE3_SCRIPTS / item.name
        if item.is_file():
            shutil.copy2(item, dest)
            print(f"    ✅ {item.name}")

print("\n✅ Reorganization complete!")
print(f"\n📁 Created directories:")
print(f"  {PHASE2_DIR}")
print(f"  {PHASE3_DIR}")
print("\n" + "="*70)
