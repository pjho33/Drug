#!/bin/bash
export PATH="/home/pjho3/miniconda3/bin:$PATH"
source /home/pjho3/miniconda3/etc/profile.d/conda.sh
cd /home/pjho3/projects/Drug
SMILES=$(cat tris_peg24_triazole_amide_l_glucose.smi)
./run_manual_pipeline.sh raw_data/4PYP_trimer.pdb "$SMILES" run_tris_peg24_triazole_amide
