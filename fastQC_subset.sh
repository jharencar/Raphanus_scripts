#!/bin/bash -l

#SBATCH --job-name=fastqc_parents_subset
#SBATCH -D /group/jrigrp11/globus-write/jgharenc/parents 
#SBATCH -e /home/jgharenc/moe_raw_seq/slurm_log/sterror_fastqc%j.txt
#SBATCH -o /home/jgharenc/moe_raw_seq/slurm_log/stdoutput_fastqc%j.txt
#SBATCH -A jrigrp
#SBATCH -p high2
#SBATCH -t 4:00:00
#SBATCH --mem 4G

module load fastqc/0.11.9

fastqc MB_P3_P2B07_S242_L003_R1_001.fastq.gz MB_P3_P2B07_S242_L003_R2_001.fastq.gz MB_P4_P1F01_S294_L003_R1_001.fastq.gz MB_P4_P1F01_S294_L003_R2_001.fastq.gz  MB_P4_P1D01_S292_L003_R1_001.fastq.gz MB_P4_P1D01_S292_L003_R2_001.fastq.gz MB_P3_P2G11_S279_L003_R2_001.fastq.gz MB_P3_P2G11_S279_L003_R1_001.fastq.gz
