#!/bin/bash
#SBATCH --partition 512x64-ib 
#SBATCH --qos=ib
#SBATCH --account=ib
#SBATCH --job-name=AHV_QTL_ALLEref_vcf_unfiltered
#SBATCH --output=AHV_QTL_ALLEref_vcf_unfiltered_.%A_%a.out
#SBATCH -e AHV_QTL_ALLEref_vcf_unfiltered_.%A_%a.err
#SBATCH -t 4-00:00:00
#SBATCH --mail-type=BEGAN,END,FAIL
#SBATCH --mail-user=jharenca@ucsc.edu
#SBATCH --mem=8G # (leaving some for operations - 4 would fully fill)
#SBATCH --cpus-per-task=30
#SBATCH --array=1-11 # set to number of scaffold.list files

# Code to generate initial vcf in Nicolas' snp filtering 3 step code; code originally written by Tyler Linderoth
# Modified by Julia Harenčár and Kate Uckele, latest update 12/14/2023

#load bcftools
module load bcftools

## change into the directory with bam.filelist (explicit paths in list) and scaffold lists
cd /hb/scratch/jharenca/QTL_mapping/VCFs

# Set up array variables 
BAM_LIST="AHV_F2_F3_bams_241121.list"
REF="/hb/groups/kay_lab/genomes/Darragh_genomes/Costus_allenii.fasta"
THREADS=60

# step 1: for each of 11 sets of scaffolds
bcftools mpileup -b ${BAM_LIST} -R scaffolds_C_alle_${SLURM_ARRAY_TASK_ID}.list -f ${REF} --threads ${THREADS} -d 100000000 -Q 0 -a AD,DP,SP,INFO/ADF,INFO/ADR -L 100000000 --skip-indels | \
bcftools call --threads ${THREADS} -m --variants-only --skip-variants indels -A -f GQ - -O z > unfiltered_ALLEref_scaffold_set${SLURM_ARRAY_TASK_ID}.vcf.gz

echo "bcftools mpileup -b ${BAM_LIST} -R scaffolds_${SLURM_ARRAY_TASK_ID}.list -f ${REF} --threads ${THREADS} -d 100000000 -Q 0 -a AD,DP,SP,INFO/ADF,INFO/ADR -L 100000000 --skip-indels | \
bcftools call --threads ${THREADS} -m --variants-only --skip-variants indels -A -f GQ - -O z > unfiltered_scaffold_set${SLURM_ARRAY_TASK_ID}.vcf.gz"

##########################################################################################
## BCFTOOLS MPILEUP OPTIONS: ##
## -d, --max-depth INT
# At a position, read maximally INT reads per input file. Note that the original samtools 
# mpileup command had a minimum value of 8000/n where n was the number of input files given 
# to mpileup. This means that in samtools mpileup the default was highly likely to be increased 
# and the -d parameter would have an effect only once above the cross-sample minimum of 8000. 
# This behavior was problematic when working with a combination of single- and multi-sample bams, 
# therefore in bcftools mpileup the user is given the full control (and responsibility), and 
# an informative message is printed instead [250]

## -Q, --min-BQ INT
# Minimum base quality for a base to be considered [13]

## -L, --max-idepth INT
# Skip INDEL calling if the average per-sample depth is above INT [250]

##########################################################################################
## BCFTOOLS CALL OPTIONS: ##
## -m, --multiallelic-caller
# alternative model for multiallelic and rare-variant calling designed to overcome known limitations in -c calling model (conflicts with -c)

## -A, --keep-alts
