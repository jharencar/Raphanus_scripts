#!/bin/bash

# === USER SETTINGS ===
ASSOC_FILE="$1"  # e.g., anthocyanins_raph1_0.2geno_gemma_output.assoc.txt
GFF_FILE="$2"    # e.g., NAU-LB.V1_annotated_fixed.gff
PVALUE_THRESHOLD=${3:-5e-08}
DISTANCE_BP=${4:-1000}

# === OUTPUT FILES ===
PEAKS_OUT="significant_peaks.tsv"
ANNOT_OUT="annotations_near_peaks.tsv"

# === STEP 1: Extract significant peaks ===
LOG_PVAL_THRESH=$(awk -v p="$PVALUE_THRESHOLD" 'BEGIN { print -log(p)/log(10) }')

awk -v thresh="$LOG_PVAL_THRESH" 'BEGIN {OFS="\t"}
    NR==1 || (-log($11)/log(10)) > thresh {
        print $1, $3, $11
    }
' "$ASSOC_FILE" > "$PEAKS_OUT"

echo "Extracted significant peaks to $PEAKS_OUT"

# === STEP 2: Convert GFF to BED ===
GFF_BED="annotations.bed"

awk 'BEGIN{FS=OFS="\t"} !/^#/ {
    match($9, /ID=([^;]+)/, arr)
    id = (arr[1] != "") ? arr[1] : "NA"
    print $1, $4 - 1, $5, id, ".", $7, $3, $9
}' "$GFF_FILE" > "$GFF_BED"

# === STEP 3: Convert peaks to BED format with flanking regions ===
PEAKS_BED="peaks.bed"
awk -v OFS="\t" -v window="$DISTANCE_BP" '
    NR > 1 {
        chr = "chr" $1
        start = $2 - window
        if (start < 0) start = 0
        end = $2 + window
        print chr, start, end, $2, $3
    }
' "$PEAKS_OUT" > "$PEAKS_BED"

# === STEP 4: Find annotations within window using bedtools ===
bedtools intersect -wa -wb -a "$PEAKS_BED" -b "$GFF_BED" > "intersections.tsv"

# === STEP 5: Format final output ===
awk 'BEGIN{FS=OFS="\t"} {
    print "peak_chr="$1, "peak_pos="$4, "pval="$5, "feature_chr="$6, "feature_start="$7+1, "feature_end="$8, \
          "feature_id="$9, "strand="$11, "type="$12, "info="$13
}' intersections.tsv > "$ANNOT_OUT"

echo "Annotations near peaks written to $ANNOT_OUT"

# === STEP 6: Extract feature sequences from reference genome ===
FASTA="/group/jrigrp11/juliagh/genomes/renamed_NAU.LB_R.sativus_genome.fasta"
FEATURES_BED="features_near_peaks.bed"
SEQUENCES_OUT="feature_sequences.fa"

# Extract only feature regions from intersected file
awk 'BEGIN{FS=OFS="\t"} {
    print $6, $7, $8, $9, ".", $11
}' intersections.tsv > "$FEATURES_BED"

# Get sequences using bedtools
bedtools getfasta -fi "$FASTA" -bed "$FEATURES_BED" -s -name > "$SEQUENCES_OUT"

echo "Feature sequences written to $SEQUENCES_OUT"
