# Phase 3: SDG-Glycated GLUT1 Complex MD Simulation

**Date**: January 2026  
**Objective**: 100ns MD simulation and MM-GBSA analysis of SDG ligand binding to glycosylated GLUT1

## Overview

This phase focuses on:
1. **Glycosylated GLUT1** - Protein with physiologically relevant glycosylation
2. **SDG Ligand** - Optimized glucose-based ligand from Phase 2
3. **Long-timescale MD** - 100ns production simulation
4. **Binding Free Energy** - MM-GBSA calculations

## Directory Structure

```
20260106_Phase3_SDG_Glycated_GLUT1/
├── scripts/              # Analysis and conversion scripts
├── md_simulation/        # MD simulation data and setup
│   ├── phase3_final/     # Final production run
│   ├── phase3_peg2_submission/
│   ├── phase3_glycosylation/
│   └── ...
├── analysis/             # Trajectory analysis results
│   ├── RMSD plots
│   ├── RMSF analysis
│   └── Hydrogen bond analysis
└── mmpbsa/              # MM-GBSA binding free energy
    ├── dry_complex.dcd   # Stripped trajectory (protein+glycan+SDG)
    ├── dry_complex.pdb
    └── analysis results
```

## System Details

- **Protein**: Glycosylated GLUT1 (with NAG, BMA, MAN, GAL, FUC, SIA glycans)
- **Ligand**: SDG (optimized from Phase 2)
- **Membrane**: POPC lipid bilayer
- **Solvent**: TIP3P water with 0.15M NaCl
- **Total atoms**: ~162,000 (full system)
- **Dry system**: ~17,000 atoms (protein + glycan + SDG)

## Simulation Protocol

1. **System Preparation**: CHARMM-GUI Membrane Builder
2. **Equilibration**: 6 steps (NVT → NPT with gradual restraint release)
3. **Production**: 100ns NPT at 310K, 1 bar
4. **Trajectory**: 10,000 frames saved

## Analysis Performed

### Trajectory Analysis
- **RMSD**: Protein and ligand stability
- **RMSF**: Per-residue flexibility
- **Hydrogen Bonds**: Protein-ligand interactions
- **Minimum Distance**: Closest contact tracking
- **COM Distance**: Center of mass separation

### Binding Free Energy (MM-GBSA)
- **Method**: Generalized Born (igb=5)
- **Salt concentration**: 0.15M
- **Frames analyzed**: 10,000
- **Components**: ΔG_bind = ΔE_MM + ΔG_sol - TΔS

## Key Scripts

### Analysis Scripts
- `analyze_trajectory_safe.py` - Memory-efficient RMSD/distance analysis
- `analyze_rmsf_hbond.py` - RMSF and hydrogen bond analysis
- `make_dry_traj.py` - Create stripped trajectory for MM-GBSA

### Conversion Scripts
- `make_dry_psf.py` - Create dry PSF (protein+glycan+ligand)
- `psf_to_amber.py` - Convert CHARMM PSF to Amber prmtop
- `add_radii.py` - Add mbondi2 radii for GB calculations

### Pipeline Scripts
- `run_mmpbsa_pipeline.sh` - Complete MM-GBSA workflow

## Data Files

### Original MD Data
Located in: `/home/pjho3tr/Downloads/GlycatedGLUT1SDGComplex260107/`
- `step5_input.psf` - Full system topology
- `step5_input.pdb` - Initial coordinates
- `production.dcd` - 100ns trajectory (~8GB)

### Processed Data
- `dry_complex.dcd` - Stripped trajectory (7,432 atoms)
- `dry_complex.pdb` - Dry system structure
- Analysis plots and CSV files

## Notes

- Glycosylation sites were identified from experimental data
- SDG ligand parameters generated with CGenFF
- All force field files from CHARMM36m
- MM-GBSA calculations use protein+glycan as receptor, SDG as ligand
