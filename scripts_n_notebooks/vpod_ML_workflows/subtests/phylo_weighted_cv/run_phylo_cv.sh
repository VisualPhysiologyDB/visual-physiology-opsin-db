#!/bin/bash  -l
#SBATCH --nodes=3 --ntasks-per-node 20
#SBATCH --time=24:00:00

conda activate deepBreaks_env

cd $SLURM_SUBMIT_DIR

python /home/seth_frazer/vpod/phylo_weighted_cv.py \
    --path ./vpod_1.2_data_splits_2024-10-19_10-30-09 \
    --seqFileName ./wt_aligned_VPOD_1.2_het.fasta \
    --metaDataFileName ./wt_meta.tsv \
    --tree_folder ./phylo_cv_results/opsin_wt_tree/vpod_1.2_wt \
    --ds wt \
    --encoding aa_prop \
    --gap_threshold 0.5 \