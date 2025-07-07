#!/bin/bash -l

#SBATCH --job-name=rm_bads_vcf
#SBATCH -D /group/jrigrp11/juliagh/angsd
#SBATCH -o /home/jgharenc/Raphanus_scripts/slurm_log/rm_bads_vcf.%A_%a.out
#SBATCH -e /home/jgharenc/Raphanus_scripts/slurm_log/rm_bads_vcf.%A_%a.err
#SBATCH -A jrigrp
#SBATCH -p high2
#SBATCH --time=10-00:00:00 #arbitratily large, no idea how long it will take
#SBATCH -c 2
#SBATCH --mem=150G

# load module
module load bcftools

bcftools view -s ^58.2,138.1 -Oz -o angsd_hardcall_250212_no58.2_138.1.vcf.gz angsd_hardcall_250214.vcf.gz
