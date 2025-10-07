#!/bin/bash  -l
#SBATCH --nodes=1 --ntasks-per-node 10
#SBATCH --time=24:00:00

cd $SLURM_SUBMIT_DIR

conda activate deepBreaks_env

python /home/seth_frazer/vpod/vpod_bootstrap_gen_one_hot_wf.py
