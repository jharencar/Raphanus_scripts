#!/bin/bash -l

#SBATCH --job-name=gemma_kin_gen
#SBATCH -D /group/jrigrp11/juliagh/GEMMA
#SBATCH -o /home/jgharenc/Raphanus_scripts/slurm_log/gemma_kin_gen.%A_%a.out
#SBATCH -e /home/jgharenc/Raphanus_scripts/slurm_log/gemma_kin_gen.%A_%a.err
#SBATCH -A jrigrp
#SBATCH -p high2
#SBATCH --time=10-00:00:00 #arbitratily large, no idea how long it will take
#SBATCH -c 1
#SBATCH --mem=500G

# load module
module load gemma/0.98.5

# genetate the kinship matrix
gemma -bfile ./input/plink_filtered_58.2_138.1_rmd_0.1geno_0.5mind -gk 1 -o raph_gemma_kinship
