#!/bin/bash -l

#SBATCH --job-name=rename_bcf
#SBATCH -D /group/jrigrp11/juliagh/angsd
#SBATCH -o /home/jgharenc/Raphanus_scripts/slurm_log/rename_bcf.%A_%a.out
#SBATCH -e /home/jgharenc/Raphanus_scripts/slurm_log/rename_bcf.%A_%a.err
#SBATCH -A jrigrp
#SBATCH -p high2
#SBATCH --time=10-00:00:00 #arbitratily large, no idea how long it will take
#SBATCH -c 2
#SBATCH --mem=150G

# load module
module load bcftools

IN_BCF=angsd_merged_bams_250214.bcf
OUT_BCF=angsd_merged_bams_renamed_250214.bcf
TMP_SAMPLE_LIST=tmp_sample_list.txt

# Extract the sample names from the BCF file
bcftools query -l "$IN_BCF" > "$TMP_SAMPLE_LIST"

# Create a new header file with the renamed samples
# This uses awk to extract the numerical part from the full path.
# Adjust the awk command if your path structure is different.
awk '{
  # Splits the line by the delimiter '/'.
  split($0, parts, "/");
  # Splits the last element of 'parts' by the delimiters '.' and '_'.
  split(parts[length(parts)], name_parts, "[._]");
  # Prints the original name and the renamed name, separated by a tab.
  print $0 "\t" name_parts[1] "." name_parts[2];
}' "$TMP_SAMPLE_LIST" > new_header.txt

# Reheader the BCF file using the new header file
bcftools reheader -s new_header.txt -o "$OUT_BCF" "$IN_BCF"

# Index the reheadered BCF file
bcftools index "$OUT_BCF"

# remove tmp files - can uncomment when sure it is working
# rm "$TMP_SAMPLE_LIST" new_header.txt

echo "Sample renaming complete.  New BCF file: $OUT_BCF"
