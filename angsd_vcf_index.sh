#!/bin/bash -l

#SBATCH --job-name=angsd_var_call
#SBATCH -D /group/jrigrp11/juliagh/moe_data_cleaned_bams
#SBATCH -o /home/jgharenc/Raphanus_scripts/slurm_log/angsd_vcf.%A_%a.out
#SBATCH -e /home/jgharenc/Raphanus_scripts/slurm_log/angsd_vcf.%A_%a.err
#SBATCH -A jrigrp
#SBATCH -p high2
#SBATCH --time=10-00:00:00 #arbitratily large, no idea how long it will take
#SBATCH -c 2
#SBATCH --mem=150G

# load module
module load bcftools

# Path to output directory
OUTDIR=/group/jrigrp11/juliagh/angsd

# Index VCF (vcf.gz) file
bcftools index $OUTDIR/*.bcf
