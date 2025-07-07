#!/bin/bash -l

#SBATCH --job-name=bcf_to_plink
#SBATCH -D /group/jrigrp11/juliagh/plink_moe_data
#SBATCH -o /home/jgharenc/Raphanus_scripts/slurm_log/bcf_to_plink.%A_%a.out
#SBATCH -e /home/jgharenc/Raphanus_scripts/slurm_log/bcf_to_plink.%A_%a.err
#SBATCH -A jrigrp
#SBATCH -p high2
#SBATCH --time=1-00:00:00 #arbitratily large, on Raph vcf.gz of 38G, took 26 min including indexing (most of that probably the indexing)
#SBATCH -c 1
#SBATCH --mem=100G

# load module
module load bcftools 

# set variables
#IN_BCF=/group/jrigrp11/juliagh/angsd/angsd_merged_bams_renamed_250214.bcf
IN_VCF=/group/jrigrp11/juliagh/angsd/angsd_hardcall_250214_no58.2_138.1.vcf.gz
MISS_IND=0.5 # set individual missingness threshold; making it lenient bc low coverage data,- e.g. a value of 0.25 excludes individuals missing data at more than 25% of sites; TOO severe! removes half of samples... removing this filter entirely for now... 
MISS_GENO=0.20 # set variant missingness threshold; e.g. 0.1 means exclude sites with more than 10% missing data
MAF=0.05 # set MAF threshold; remove sites with minor allele freq less than $MAF
HWE=1e-6 # Set HWE threshold; filters SNPs that deviate significantly from Hardy-Weinberg equilibrium (HWE) at a significance level of 1e-6 (maybe adjust? not sure how to select this sig level).

# NOTE: In PLINK, the HWE calculation is done *before* the MAF filtering. To filter based on HWE more stringently, do so separately after MAF filtering.

# Convert BCF to VCF (PLINK can read VCF but BCFtools is often faster)
#bcftools view "$IN_BCF" -o "$IN_VCF" -O z

# Index the VCF
#bcftools index "$IN_VCF"

module load conda/plink/1.90b6.21

# Convert VCF to PLINK 1.9 format, filter out indels, Missingness filter (SNPs and individuals), MAF filter, and HWE filter
plink --vcf "$IN_VCF" \
      --snps-only just-acgt \
      --geno $MISS_GENO \
      --maf $MAF \
      --hwe $HWE \
      --keep-allele-order \
      --allow-extra-chr \
      --make-bed \
      --out temp_plink_unfiltered 

# filter individuals with less than 50% of the filtered sites from above covered

# Missingness filter (SNPs and individuals), MAF filter, and HWE filter
plink --bfile temp_plink_unfiltered \
      --mind $MISS_IND \
      --make-bed \
      --keep-allele-order \
      --allow-extra-chr \
      --out plink_filtered_58.2_138.1_rmd_0.2geno_0.5mind 
