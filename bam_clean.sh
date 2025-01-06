#!/bin/bash -l

#SBATCH -D /group/jrigrp11/juliagh/moe_data_deduped_bams
#SBATCH -o /home/jgharenc/Raphanus_scripts/slurm_log/clean_bams.%A_%a.txt
#SBATCH -e /home/jgharenc/Raphanus_scripts/slurm_log/clean_bams.%A_%a.txt
#SBATCH -A jrigrp
#SBATCH -p high2
#SBATCH --time=10-00:00:00
#SBATCH --array=1-382 # []%6 limits it so only 6 run in parallel at once
#SBATCH --mem=60G # Maybe bump up to like 20 and see if it runs faster?

# Bam cleaning code: removes duplicates, low mapping quality reads, and regions with greater or less than a 99th percentile coverage cutoff
#by: Julia Harenčár, 210817

## load samtools
module load samtools/1.16.1

# set variables
BAM_FILE=$(ls *_paired_output.bam | sed -n ${SLURM_ARRAY_TASK_ID}p) # make list of file names to iterate through
## FIX BELOW BASED ON FILE NAMES ##
SAMPLE_ID=$(echo $BAM_FILE |cut  -d "_" -f 1,2)

############################################################
### STEP 1: REMOVING DUPLICATES and FIXING MATE INFO ###
# to use markdup, need several steps first:

#name order/collate
#echo "samtools collate -o ${SAMPLE_ID}_collated.bam ${BAM_FILE}" 
#samtools collate -o ${SAMPLE_ID}_collated.bam ${BAM_FILE} || echo "samtools collate failed for ${SAMPLE_ID}"

#Add ms and MC tags for markdup to use later (-m adds ms (mate score) tags):
#samtools fixmate -m ${SAMPLE_ID}_collated.bam ${SAMPLE_ID}_fixmated.bam || echo "samtools fixmate failed for ${SAMPLE_ID}"

#Markdup needs position order (-o is output file)
#samtools sort -o ${SAMPLE_ID}_sorted.bam ${SAMPLE_ID}_fixmated.bam || echo "samtools sort1 failed for ${SAMPLE_ID}"

#finally, mark and remove duplicates
#samtools markdup -r ${SAMPLE_ID}_sorted.bam ${SAMPLE_ID}_rmdup.bam || echo "samtools markdup failed for ${SAMPLE_ID}"

#remove all the annoying intermediate bams
#rm ${SAMPLE_ID}_collated.bam ${SAMPLE_ID}_sorted.bam ${SAMPLE_ID}_fixmated.bam 

############################################################
### STEP 2: REMOVING LOW MAP QUALITY READS ###

# indexing/sorting bams
#samtools index ${SAMPLE_ID}_rmdup.bam || echo "samtools index1 failed for ${SAMPLE_ID}"

# Filter bams to remove reads with a MAPQ phred score less than 20
#samtools view -h -q 20 ${SAMPLE_ID}_rmdup.bam > ${SAMPLE_ID}_rmdup_q20.bam || echo "samtools q20 filter failed for ${SAMPLE_ID}"

############################################################
### STEP 3: REMOVING READS IN 99TH PERCENTILE OF COVERAGE

# calculate coverage
#samtools depth -a ${SAMPLE_ID}_rmdup_q20.bam > ${SAMPLE_ID}_coverage.txt || echo "samtools depth failed for ${SAMPLE_ID}"

# load conda env with python dependencies
module load conda
conda activate datasci

# run python section for calculating threshold cutoffs using a heredoc
python3 <<EOF
# get 99th percentile thresholds with python
import numpy as np
import pandas as pd

# Define SAMPLE_ID from Bash variable
SAMPLE_ID = "${SAMPLE_ID}"

# Load coverage file
coverage = pd.read_csv("${SAMPLE_ID}_coverage.txt", sep="\t", header=None, names=["chrom", "pos", "cov"])

# Calculate thresholds
low_thresh = np.percentile(coverage['cov'], 1)  # 1st percentile
high_thresh = np.percentile(coverage['cov'], 99)  # 99th percentile

print(f"Low coverage threshold: {low_thresh}")
print(f"High coverage threshold: {high_thresh}")

# Filter for positions outside thresholds
exclude = coverage[(coverage['cov'] < low_thresh) | (coverage['cov'] > high_thresh)].copy()

# Group consecutive positions into regions (continuous stretches)
exclude['start'] = exclude['pos'] - 1  # BED format is 0-based, so convert 'start' to 0-based 'pos'
exclude['gap'] = (exclude['pos'] != exclude['pos'].shift() + 1).astype(int).cumsum() # if the current position has a gap >1 with the previous condition, result = true meaning there is a gap, converts t/f to 0/1, then cumsum is the 'group number', where consecutive pos are all in a group, group changes after a gap
bed = exclude.groupby(['chrom', 'gap']).agg({'start': 'min', 'pos': 'max'}).reset_index() #Combines positions in same group (no gap) into a single region

# Save to BED file
bed[['chrom', 'start', 'pos']].rename(columns={'pos': 'end'}).to_csv(f"{SAMPLE_ID}_exclude.bed", sep="\t", header=False, index=False)
EOF

# Filter the bam based on the BED file with regions to exclude
samtools view -h -L ${SAMPLE_ID}_exclude.bed -U ${SAMPLE_ID}_rmdup_q20_99covCO.bam ${SAMPLE_ID}_rmdup_q20.bam > /dev/null || echo "samtools filter by BED failed for ${SAMPLE_ID}"

# Sort final cleaned bams
samtools sort -o ${SAMPLE_ID}_rmdup_q20_99covCO_sorted.bam ${SAMPLE_ID}_rmdup_q20_99covCO.bam || echo "samtools sort2 failed for ${SAMPLE_ID}"

# Index final cleaned bams
samtools index ${SAMPLE_ID}_rmdup_q20_99covCO_sorted.bam || echo "samtools index2 failed for ${SAMPLE_ID}"
