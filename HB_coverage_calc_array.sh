#!/bin/bash
#SBATCH --partition=512x64-ib 
#SBATCH --qos=ib     
#SBATCH --account=ib 
#SBATCH --mail-type=ALL              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jharenca@ucsc.edu  # Where to send mail
#SBATCH --job-name=bam_coverage # Job name
#SBATCH --output=bam_coverage_.%A_%a.out
#SBATCH -e bam_coverage_.%A_%a.err
#SBATCH --ntasks-per-node=10
#SBATCH --array=1-442  # []%3 limits it so only 3 run in parallel at once
#submit with SBATCH array, submit as normal or to specify certain array numbers: sbatch --array=x-y job_script.sbatch

## load samtools
module load samtools

## change into the directory with .bams
cd /hb/scratch/jharenca/QTL_mapping/cleaning/trimmed

# variable setting
BAM_FILE=$(ls *_paired_output.bam | sed -n ${SLURM_ARRAY_TASK_ID}p) # make list of file names to iterate through

# EG: name pulls 19.572 and SPP pulls LAEV from: 19.572_LAEV_merged.bam
SAMPLE_ID=$(echo $BAM_FILE |cut  -d "_" -f 1,2 )
#SAMPLE_SPP=$(echo $BAM_FILE |cut  -d "_" -f 2 )

# Get mean read depth - here 's' is cumulative depth at all sites and 'c' is the total # of positions
echo " samtools depth -a ${BAM_FILE} | awk '{c++;s+=$3}END{print s/c}' "
depth=$(samtools depth -a ${BAM_FILE} | awk '{c++;s+=$3}END{print s/c}') 

# Get Breadth of Coverage - here 'c' is the total # of positions again and total receives an increment of 1 when depth > 0
# bredth is then (total/c)*100 (percent of the genome covered at least once!)
echo " samtools depth -a ${BAM_FILE} | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' "
breadth=$(samtools depth -a ${BAM_FILE} | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}')

# Print depth and breadth with file name to .csv
echo "${SAMPLE_ID}, $depth, $breadth" >> qtl_mapping_pop_coverage_data_PREclean.csv
