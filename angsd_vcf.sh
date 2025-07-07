#!/bin/bash -l

#SBATCH --job-name=angsd_var_call
#SBATCH -D /group/jrigrp11/juliagh/moe_data_cleaned_bams/bams_merged_by_id
#SBATCH -o /home/jgharenc/Raphanus_scripts/slurm_log/angsd_vcf.%A_%a.out
#SBATCH -e /home/jgharenc/Raphanus_scripts/slurm_log/angsd_vcf.%A_%a.err
#SBATCH -A jrigrp
#SBATCH -p high2
#SBATCH --time=10-00:00:00 #arbitratily large, no idea how long it will take
#SBATCH -c 2
#SBATCH --mem=150G

# load module
module load angsd/0.935 bcftools

# set variables
REF="/group/jrigrp11/juliagh/genomes/renamed_NAU.LB_R.sativus_genome.fasta.gz"
BAMLIST="clean_merged_bams_250214.list"

# Path to output directory
OUTDIR=/group/jrigrp11/juliagh/angsd

# Call SNPs and genotypes with ANGSD
angsd \
-nThreads 2 \
-bam $BAMLIST \
-ref $REF \
-out $OUTDIR/angsd_merged_bams_250214 \
-minMapQ 20 \
-minQ 20 \
-dobcf 1 \
-gl 1 \
-doMaf 1 \
-doMajorMinor 1 \
-SNP_pval 1e-6 \
-doPost 1 \
-doCounts 1 \
-doGeno 1

# Index VCF (vcf.gz) file
bcftools index $OUTDIR/*.vcf.gz
