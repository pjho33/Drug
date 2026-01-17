# Drug Discovery Project: GLUT1 Inhibitors

**Project Goal**: Design and validate glucose-based inhibitors for GLUT1 transporter

## Project Timeline

### Phase 2: Initial Ligand Screening (December 2024)
**Folder**: `20241209_Phase2_Glucose_BNG_Tripod/`

Tested multiple glucose-based scaffolds:
- D-glucose and L-glucose
- BNG derivatives
- Tripodal scaffolds with varying PEG linkers (PEG2, PEG6, PEG12, PEG24)
- Tris-based scaffolds

**Key Findings**:
- 3 replicas per ligand for reproducibility
- Identified optimal PEG linker lengths
- Selected SDG as lead compound

### Phase 3: Lead Optimization (January 2026)
**Folder**: `20260106_Phase3_SDG_Glycated_GLUT1/`

Deep characterization of SDG binding to glycosylated GLUT1:
- 100ns MD simulation
- Glycosylated GLUT1 (physiologically relevant)
- MM-GBSA binding free energy calculations
- Detailed interaction analysis

**Key Results**:
- Stable binding throughout 100ns
- Favorable binding free energy
- Key hydrogen bonds identified
- Glycan interactions characterized

## Directory Structure

```
Drug/
├── 20241209_Phase2_Glucose_BNG_Tripod/
│   ├── scripts/          # Docking and MD pipeline
│   ├── results/          # All replica simulations
│   ├── analysis/         # Analysis results
│   └── structures/       # Ligand structures
│
├── 20260106_Phase3_SDG_Glycated_GLUT1/
│   ├── scripts/          # Analysis and conversion scripts
│   ├── md_simulation/    # 100ns MD data
│   ├── analysis/         # Trajectory analysis
│   └── mmpbsa/          # Binding free energy
│
├── scripts/              # Current working scripts
└── mmpbsa/              # Current MM-GBSA analysis
```

## Methodology

### Molecular Docking
- **Tools**: GNINA, DiffDock
- **Receptor**: GLUT1 (PDB-based, glycosylated)
- **Validation**: Redocking known inhibitors

### Molecular Dynamics
- **Software**: OpenMM 8.0
- **Force Field**: CHARMM36m
- **System Setup**: CHARMM-GUI Membrane Builder
- **Membrane**: POPC lipid bilayer
- **Solvent**: TIP3P water, 0.15M NaCl
- **Temperature**: 310K (physiological)
- **Pressure**: 1 bar

### Binding Free Energy
- **Method**: MM-GBSA (Molecular Mechanics - Generalized Born Surface Area)
- **Tool**: MMPBSA.py (AmberTools)
- **GB Model**: igb=5 (mbondi2 radii)
- **Analysis**: Per-frame decomposition

## Key Technologies

- **CHARMM-GUI**: System preparation
- **OpenMM**: MD simulations
- **MDAnalysis**: Trajectory analysis
- **ParmEd**: Topology conversion (CHARMM ↔ Amber)
- **AmberTools**: MM-GBSA calculations
- **Python**: Analysis automation

## Data Management

### Original Data
- Phase 2: `/home/pjho3tr/projects/Drug_local_20160117/`
- Phase 3 MD: `/home/pjho3tr/Downloads/GlycatedGLUT1SDGComplex260107/`

### Organized Data
- All results reorganized into dated folders
- README files in each folder
- Scripts preserved with results

## Future Directions

1. **Control Simulations**: Non-glycosylated GLUT1
2. **Comparative Analysis**: Glycosylated vs. non-glycosylated
3. **Experimental Validation**: Binding assays
4. **Optimization**: Further SDG modifications

## References

- CHARMM-GUI: https://www.charmm-gui.org/
- OpenMM: https://openmm.org/
- MDAnalysis: https://www.mdanalysis.org/
- AmberTools: https://ambermd.org/AmberTools.php

---

**Last Updated**: January 17, 2026
