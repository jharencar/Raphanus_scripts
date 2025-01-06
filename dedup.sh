#!/bin/bash -l

#SBATCH --job-name=dedup # Job name
#SBATCH -D /group/jrigrp11/globus-write/jgharenc/tmp_undeduped/ 
#SBATCH -e /home/jgharenc/Raphanus_scripts/slurm_log/dedup%j.txt
#SBATCH -o /home/jgharenc/Raphanus_scripts/slurm_log/dedup%j.txt
#SBATCH -A jrigrp
#SBATCH -p high2
#SBATCH --time=1-12:00:00  
#SBATCH --mem=80G                    # Memory per node - was at 20G when first run on farm on 12/2/24; may need more based on initial output

# load modules
module load openjdk

## call the code to run
python3 /home/jgharenc/Raphanus_scripts/dedup.py
