#!/bin/bash -l

#SBATCH --job-name=merge_bams_by_ID
#SBATCH -D /group/jrigrp11/juliagh/moe_data_cleaned_bams
#SBATCH -o /home/jgharenc/Raphanus_scripts/slurm_log/merge_bams.%A_%a.out
#SBATCH -e /home/jgharenc/Raphanus_scripts/slurm_log/merge_bams.%A_%a.err
#SBATCH -A jrigrp
#SBATCH -p high2
#SBATCH --time=10-00:00:00 #arbitratily large, no idea how long it will take
#SBATCH -c 8 # selecting 4 because of samtools merge -@ 8, but not sure that will work that simply... (might need 9 bc i think it is 4 in addition to the main?)
#SBATCH --mem=150G # just a random guess...

# load modules and acivate conda env with pandas
module load samtools conda 
conda activate my_env

# run python script to merge
python /home/jgharenc/Raphanus_scripts/merge_bams_by_id.py 
