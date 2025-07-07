#!/bin/bash
#SBATCH --job-name=extract_GEMMA_peaks_annotations
#SBATCH -D /group/jrigrp11/juliagh/GEMMA/output/log_response_ratio_0.2geno 
#SBATCH -e /home/jgharenc/Raphanus_scripts/slurm_log/extract_peaks_annotations.%j.txt
#SBATCH -o /home/jgharenc/Raphanus_scripts/slurm_log/extract_peaks_annotations.%j.txt
#SBATCH -A jrigrp
#SBATCH -p high2
#SBATCH --time=20:00:00
#SBATCH --mem=240Gb

# load bedtools
module load bedtools2/2.31.1   

# run the script specifying input files and peak threshold and annotation window size
/home/jgharenc/Raphanus_scripts/extract_gemma_peaks_and_nearby_annotations.sh anthocyanins_0.2geno_gemma_output.assoc.txt NAU-LB.V1_annotated_fixed.gff 5e-08 3000
