#!/bin/bash
#SBATCH --job-name=bwa_mapping_R.sativus 
#SBATCH -D /group/jrigrp11/juliagh/moe_data_deduped_bams/tmp/ 
#SBATCH -o /home/jgharenc/Raphanus_scripts/slurm_log/bwa_mapping.%A_%a.txt
#SBATCH -e /home/jgharenc/Raphanus_scripts/slurm_log/bwa_mapping.%A_%a.txt
#SBATCH -c 40
#SBATCH -A jrigrp
#SBATCH -p high2
#SBATCH --time=1-00:00:00
#SBATCH --array=1-383 # []%6 limits it so only 6 run in parallel at once
#submit with SBATCH array, submit as normal or to specify certain array numbers: sbatch --array=x-y job_script.sbatch

##(more info on arrays: https://support.ceci-hpc.be/doc/_contents/SubmittingJobs/SlurmFAQ.html)

##get bwa mem and samtools modules
module load bwa-mem2 samtools

# variable setting
NUMTHREADS=40
LIBRARY_NAME="moe_seqs"
REF="/group/jrigrp11/juliagh/genomes/NAU.LB_R.sativus_genome.fasta.gz"

# Set up array variables so that the same SLURM_ARRAY_TASK_ID references both input files for a given individual (F and R reads)
R1=$(ls *_clumped1.fq.gz| sed -n ${SLURM_ARRAY_TASK_ID}p)
R2=$(ls *_clumped2.fq.gz| sed -n ${SLURM_ARRAY_TASK_ID}p)

# pull sample name from file name
SAMPLE_NAME=$(echo $R1 |cut  -d "_" -f 2,3,4 )

### NOTE: index reference before submitting SLURM: bwa index $REF 
# bwa mapping; -t number of treads; -M is a verboseLevel ; -R readGroupHeader ; 
echo "bwa-mem2 mem -t ${NUMTHREADS} -M -R "@RG\tID:${LIBRARY_NAME}\tSM:${SAMPLE_NAME}\tLB:${LIBRARY_NAME}\tPL:ILLUMINA\tPU:none" ${REF} ${R1} ${R2} | samtools view -hb -@ ${NUMTHREADS} - | samtools sort -@ ${NUMTHREADS} - -o ${SAMPLE_NAME}_paired_output.bam"
bwa-mem2 mem -t ${NUMTHREADS} -M -R "@RG\tID:${LIBRARY_NAME}\tSM:${SAMPLE_NAME}\tLB:${LIBRARY_NAME}\tPL:ILLUMINA\tPU:none" ${REF} ${R1} ${R2} | samtools view -hb -@ ${NUMTHREADS} - | samtools sort -@ ${NUMTHREADS} - -o ${SAMPLE_NAME}_sativus_ref_paired_output.bam
