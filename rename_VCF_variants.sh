#!/bin/bash -l

#SBATCH -D /group/jrigrp11/juliagh/
#SBATCH -o /home/jgharenc/Raphanus_scripts/slurm_log/variant_name.%A_%a.out
#SBATCH -e /home/jgharenc/Raphanus_scripts/slurm_log/variant_name.%A_%a.err
#SBATCH -A jrigrp
#SBATCH -p high2
#SBATCH --time=5:00:00 #arbitratily large, no idea how long it will take
#SBATCH -c 1
#SBATCH --mem=5G

module load bcftools

#unzip vcf:
gunzip /group/jrigrp11/juliagh/angsd/angsd_hardcall_250214.vcf.gz

# set variables
IN_VCF=/group/jrigrp11/juliagh/angsd/angsd_hardcall_250214.vcf
OUT_VCF=/group/jrigrp11/juliagh/angsd/angsd_hardcall_250214_test.vcf

awk -F'\t' -v OFS='\t' '{if (!/^#/ && $3 == ".") $3 = $1"_"$2; print}' "$IN_VCF" > "$OUT_VCF" 

bgzip "$OUT_VCF"
