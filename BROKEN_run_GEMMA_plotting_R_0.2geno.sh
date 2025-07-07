#!/bin/bash
#SBATCH --job-name=plot_GEMMA
#SBATCH -D /group/jrigrp11/juliagh/GEMMA/output/log_response_ratio_0.2geno 
#SBATCH -e /home/jgharenc/Raphanus_scripts/slurm_log/plot_GEMMA%j.txt
#SBATCH -o /home/jgharenc/Raphanus_scripts/slurm_log/plot_GEMMA%j.txt
#SBATCH -A jrigrp
#SBATCH -p high2
#SBATCH --time=20:00:00
#SBATCH --mem=240Gb

module load R

# run the code - remember to first open the script and change paths to make sure it is making figures from the right files
Rscript /home/jgharenc/Raphanus_scripts/plot_GEMMA_output_0.2geno.R > $SLURM_JOBID.out
