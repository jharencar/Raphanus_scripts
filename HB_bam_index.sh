#!/bin/bash
#SBATCH --partition=128x24
#SBATCH --job-name=index # Job name
#SBATCH --time=72:00:00
#SBATCH --output=index.out
#SBATCH -e index.err
#SBATCH --mail-type=END,FAIL              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jharenca@ucsc.edu  # Where to send mail
#SBATCH --mem=10G

#load samtools
module load samtools

## change into the directory with .bams
cd /hb/scratch/jharenca/QTL_mapping/cleaning/trimmed

# set variables
BAM_LIST=$(ls *q20.bam) # make list of file names to iterate through

for BAM_NAME in $BAM_LIST; do
	samtools index ${BAM_NAME}
done
