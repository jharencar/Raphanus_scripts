#!/bin/bash -l

#SBATCH --job-name=gemma_kinship_decomp
#SBATCH -D /group/jrigrp11/juliagh/GEMMA
#SBATCH -o /home/jgharenc/Raphanus_scripts/slurm_log/gemma_kinship_decomp.%A_%a.out
#SBATCH -e /home/jgharenc/Raphanus_scripts/slurm_log/gemma_kinship_decomp.%A_%a.err
#SBATCH -A jrigrp
#SBATCH -p high2
#SBATCH --time=10-00:00:00 #arbitratily large, no idea how long it will take
#SBATCH -c 1
#SBATCH --mem=150G

# load module
module load gemma/0.98.5

# do eigen-decomposition of the relatedness matrix - help get an idea of how big a factor the relatedness is
#gemma -bfile ./input/plink_filtered_58.2_138.1_rmd -k ./output/raph_gemma_kinship.cXX.txt -eigen -o kinship_matrix_eigen-decomp

# see how much of variance in anthos explained by kinship:
gemma -bfile ./input/plink_filtered_58.2_138.1_rmd -k ./output/raph_gemma_kinship.cXX.txt -lmm 1 -o null_model_antho
