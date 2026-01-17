# Phase 2: Glucose, BNG, Tripod MD Simulations

**Date**: December 2024  
**Objective**: Molecular dynamics simulations and analysis of glucose-based ligands binding to GLUT1

## Overview

This phase focused on testing three types of ligands:
- **Glucose** (D-glucose, L-glucose)
- **BNG** (Glucose derivatives)
- **Tripod** (Tripodal glucose scaffolds with various PEG linkers)

Each ligand was tested with **3 replicas** to ensure reproducibility.

## Directory Structure

```
20241209_Phase2_Glucose_BNG_Tripod/
├── scripts/          # Pipeline scripts for docking and MD
├── results/          # MD simulation results for all ligands
│   ├── phase2_rep1/  # Replica 1 results
│   ├── phase2_rep2/  # Replica 2 results
│   ├── phase2_rep3/  # Replica 3 results
│   └── ...           # Other test runs
├── analysis/         # Analysis scripts and results
└── structures/       # Ligand structures (SMILES, SDF, PDB)
```

## Ligand Variants Tested

### Tripod Series
- `tripod_l_glucose.smi` - Basic tripod with L-glucose
- `tripod_peg2_l_glucose.smi` - Tripod with PEG2 linkers
- `tripod_peg6_l_glucose.smi` - Tripod with PEG6 linkers
- `tripod_peg12_l_glucose.smi` - Tripod with PEG12 linkers
- `tripod_peg24_l_glucose.smi` - Tripod with PEG24 linkers

### Tris Series
- `tris_peg12_l_glucose.smi` - Tris scaffold with PEG12
- `tris_peg24_triazole_amide_l_glucose.smi` - Tris with triazole-amide linker

## Key Results

Results are organized by replica and ligand type in the `results/` folder.

## Pipeline Scripts

- `run_pipeline.sh` - Main automated pipeline
- `run_manual_pipeline.sh` - Manual step-by-step pipeline
- `run_tris_peg12.sh` - Specific script for tris-PEG12 variant

## Notes

- All simulations used CHARMM-GUI for system preparation
- MD simulations performed with OpenMM
- Analysis includes RMSD, RMSF, hydrogen bonds, and binding energy calculations
