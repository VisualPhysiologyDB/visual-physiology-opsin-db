#!/bin/bash  -l
#SBATCH --nodes=3 --ntasks-per-node 10
#SBATCH --time=36:00:00

cd $SLURM_SUBMIT_DIR

conda activate deepBreaks_env

python /home/seth_frazer/vpod/vpod_bootstrap_gen_aa_prop_wf.py
