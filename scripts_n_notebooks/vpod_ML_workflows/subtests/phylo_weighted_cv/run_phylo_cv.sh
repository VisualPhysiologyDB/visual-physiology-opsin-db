#!/bin/bash  -l
#SBATCH --nodes=3 --ntasks-per-node 10
#SBATCH --time=36:00:00

conda activate deepBreaks_env

cd $SLURM_SUBMIT_DIR

python /home/seth_frazer/vpod/phylo_weighted_cv.py \
    --path ./vpod_1.2_data_splits_2025-02-28_15-51-04 \
    --seqFileName wt_aligned_VPOD_1.2_het.fasta \
    --metaDataFileName wt_meta.tsv \
    --tree_folder /home/seth_frazer/vpod/phylo_cv/trees/vpod_1.2/vpod_1.2_wt \
    --ds wt \
    --encoding aa_prop \