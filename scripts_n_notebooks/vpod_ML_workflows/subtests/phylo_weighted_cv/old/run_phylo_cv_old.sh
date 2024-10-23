#!/bin/bash  -l
#SBATCH --nodes=3 --ntasks-per-node 20
#SBATCH --time=24:00:00

conda activate deepBreaks_env

cd $SLURM_SUBMIT_DIR

python /home/seth_frazer/vpod/phylo_weighted_cv_iterate_aa_prop2.py
