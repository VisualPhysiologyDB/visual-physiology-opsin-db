#!/bin/bash  -l
#SBATCH --nodes=3 --ntasks-per-node 15
#SBATCH --time=100:00:00

cd $SLURM_SUBMIT_DIR

conda activate deepBreaks_env

python /home/seth_frazer/vpod/vpod_aa_prop_grid_search_iter.py
