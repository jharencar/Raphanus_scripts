#!/bin/bash -l

#SBATCH --job-name=merge_bams_by_ID
#SBATCH -D /group/jrigrp11/juliagh/moe_data_cleaned_bams
#SBATCH -o /home/jgharenc/Raphanus_scripts/slurm_log/merge_bams.%A_%a.out
#SBATCH -e /home/jgharenc/Raphanus_scripts/slurm_log/merge_bams.%A_%a.err
#SBATCH -A jrigrp
#SBATCH -p high2
#SBATCH --time=10-00:00:00 #arbitratily large, no idea how long it will take
#SBATCH -c 4 # selecting 4 because of samtools merge -@ 4, but not sure that will work that simply... (might need 5 bc i think it is 4 in addition to the main?)
#SBATCH --mem=150G # just a random guess...

# Exit immediately if a command exits with a non-zero status (helps debugging)
set -e

# load module
module load samtools

# Define input CSV and BAM directory
CSV_FILE="samp_to_seq_ID.csv"  
BAM_DIR="/group/jrigrp11/juliagh/moe_data_cleaned_bams" # Directory containing BAM files
OUTPUT_DIR="/group/jrigrp11/juliagh/moe_data_merged_by_indiv_bams"    # Directory for merged BAMs

# echos to ensure things are right
echo "CSV_FILE: $CSV_FILE"
echo "BAM_DIR: $BAM_DIR"
echo "OUTPUT_DIR: $OUTPUT_DIR"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Declare an associative array to store BAM files by sample_ID
declare -A sample_bams

# Read the CSV file (skipping the header with -n +2)
shopt -s nullglob #Important to set to prevent samtools merge from failing due to failed glob patterns
tail -n +2 "$CSV_FILE" | awk -F ',' '{print $1, $2}' | while read sample_id seq_id; do
	echo "Processing sample_id: $sample_id, seq_id: $seq_id"

    # Find the corresponding BAM file based on seq_ID (2>/dev/null suppresses an error if the bam is missing as that is checked later and will print a warning)
   bam_file=$(find "$BAM_DIR" -maxdepth 1 -name "${seq_id}"*_rmdup_q20_99covCO_sorted_chrs_only.bam)
 
    # Check if the BAM file exists
    if [[ -n "$bam_file" ]]; then
        sample_bams[$sample_id]="${sample_bams[$sample_id]} $bam_file" # Append bams in the associative array 'sample_bams'
        echo "Sample $sample_id -> BAM(s): $bam_file" 
    else
        echo "Warning: No BAM file found for seq_ID $seq_id (Sample $sample_id)"
    fi
done
shopt -u nullglob #Unset nullglob #reset this shell option

echo "Finished reading CSV file."
# Print the contents of the associative array for debugging
echo "Contents of sample_bams array:"
for sample in "${!sample_bams[@]}"; do
  echo "Sample ID: $sample, BAM files: ${sample_bams[$sample]}"
done

# Merge BAMs by sample_ID (-@ 4 indicates number of compression threads to use in addition to main thread)
for sample in "${!sample_bams[@]}"; do
    output_bam="$OUTPUT_DIR/${sample}_moe_ID.bam"
    echo "Merging BAMs for sample $sample..."
    
    #Samtools merge does not like leading or trailing spaces; xargs removes that whitespace 
    bams_to_merge=$(echo -n "${sample_bams[$sample]}" | xargs)
    
    samtools merge -@ 4 "$output_bam" ${sample_bams[$sample]}
    # Check if the merge was successful 
    if [ $? -eq 0 ]; then
	    echo "Successfully merged BAMs for sample $sample" 
	    echo "Indexing $output_bam" 
	    samtools index "$output_bam" 
	else 
		echo "ERROR: samtools merge failed for sample $sample" >&2 
	fi 
done

echo "Merging completed!"
