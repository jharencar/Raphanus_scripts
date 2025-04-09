#!/bin/bash -l  

#SBATCH --job-name=rename_chrs                             
#SBATCH -D /group/jrigrp11/juliagh/moe_data_deduped_bams/tmp 
#SBATCH -e /home/jgharenc/Raphanus_scripts/slurm_log/rehead_bams.%A_%a.err
#SBATCH -o /home/jgharenc/Raphanus_scripts/slurm_log/rehead_bams.%A_%a.out
#SBATCH -A jrigrp                                                     
#SBATCH -p high2                                                      
#SBATCH -t 5-00:00:00                                                    
#SBATCH --mem=10G
#SBATCH --array=1-2

# load samtools
module load samtools/1.16.1 

# set variables
BAM_FILE=$(ls *99covCO_sorted.bam | sed -n ${SLURM_ARRAY_TASK_ID}p)
SAMPLE_ID=$(echo $BAM_FILE | cut -d "_" -f 1,2)

# extract current header
echo "samtools view -H ${BAM_FILE} > ${SAMPLE_ID}_original_header.sam"
samtools view -H ${BAM_FILE} > ${SAMPLE_ID}_original_header.sam

# change chromosome names to make new header file
echo "sed -E 's/GWHCBIT0000000([1-9])/chr\1/g' ${SAMPLE_ID}_original_header.sam > ${SAMPLE_ID}_new_header.sam"
sed -E 's/GWHCBIT0000000([1-9])/chr\1/g' ${SAMPLE_ID}_original_header.sam > ${SAMPLE_ID}_new_header.sam

# rehead the bam
echo "samtools reheader ${SAMPLE_ID}_new_header.sam ${BAM_FILE} > ${SAMPLE_ID}_rmdup_q20_99covCO_sorted_renamed.bam"
samtools reheader ${SAMPLE_ID}_new_header.sam ${BAM_FILE} > ${SAMPLE_ID}_rmdup_q20_99covCO_sorted_renamed.bam

# index new bam for next step
echo "samtools index ${SAMPLE_ID}_rmdup_q20_99covCO_sorted_renamed.bam"
samtools index ${SAMPLE_ID}_rmdup_q20_99covCO_sorted_renamed.bam

# remove non-chromosome scaffolds
echo "samtools view -b ${SAMPLE_ID}_rmdup_q20_99covCO_sorted_renamed.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 > ${SAMPLE_ID}_rmdup_q20_99covCO_sorted_chrs_only.bam"
samtools view -b ${SAMPLE_ID}_rmdup_q20_99covCO_sorted_renamed.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 > ${SAMPLE_ID}_rmdup_q20_99covCO_sorted_chrs_only.bam

# index the new bams
echo "samtools index ${SAMPLE_ID}_rmdup_q20_99covCO_sorted_chrs_only.bam"
samtools index ${SAMPLE_ID}_rmdup_q20_99covCO_sorted_chrs_only.bam
