#!/bin/bash -l  

#SBATCH --job-name=qualimap_subset                             
#SBATCH -D /group/jrigrp11/juliagh/moe_data_deduped_bams/qualimap_subset 
#SBATCH -e /home/jgharenc/Raphanus_scripts/slurm_log/sterror_qualimap%j.txt  
#SBATCH -o /home/jgharenc/Raphanus_scripts/slurm_log/stdoutput_qualimap%j.txt
#SBATCH -A jrigrp                                                     
#SBATCH -p high2                                                      
#SBATCH -t 24:00:00                                                    
#SBATCH --mem=105G
#SBATCH --array=1-5

#load module
module load qualimap openjdk

# set variables
BAM_FILE=$(ls *q20.bam | sed -n ${SLURM_ARRAY_TASK_ID}p) # make list of file names to iterate through
SAMPLE_ID=$(echo $BAM_FILE |cut  -d "_" -f 1,2)

# run qualimap
qualimap bamqc -bam ${BAM_FILE} \
	--java-mem-size=100G \
	-outfile ${SAMPLE_ID}_qualimap.pdf \
	-outformat PDF

#-sd means don't count duplicates when computing coverage
#-sdmode 0 means when skipping duplicates, only skip reads flagged as duplicates
#-os means calculate stats for regions outside the feature file (I gave it a bed file "new.bed" containing coordinates of coding exons. It will then report coverage statistics for both inside and outside of coding regions)
#--java-mem-size=4500M I think I had issues of it running too slowly. You can speed it up by specifying the amount of memory available. For this script, I specified 50GB of memory in my hummingbird bash script, so I specified 45000M (50GB) for bamqc to run (I think you're supposed to do a little less than the total amount you have available).
