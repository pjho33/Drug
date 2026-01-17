#!/bin/bash
source ~/.bashrc
cd /home/pjho3/projects/Drug

SMILES=$(cat tris_peg12_l_glucose.smi | head -1)
bash run_manual_pipeline.sh raw_data/4PYP_trimer.pdb "$SMILES" run_tris_peg12
