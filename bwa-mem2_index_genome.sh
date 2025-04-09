#!/bin/bash
#SBATCH --job-name=index
#SBATCH -D /group/jrigrp11/juliagh/genomes/ 
#SBATCH -e /home/jgharenc/Raphanus_scripts/slurm_log/index%j.txt
#SBATCH -o /home/jgharenc/Raphanus_scripts/slurm_log/index%j.txt
#SBATCH -A jrigrp
#SBATCH -p high2
#SBATCH --time=04:00:00
#SBATCH --mem=240Gb

module load bwa-mem2/2.2.1 

/home/jgharenc/software/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index renamed_NAU.LB_R.sativus_genome.fasta.gz
