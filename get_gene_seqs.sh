#!/bin/bash

# load bedtools
module load bedtools2/2.31.1

# === USER INPUT ===
GFF_FILE="$1"       # e.g., chr1_6074654_6155392.gff
FASTA_FILE="/group/jrigrp11/juliagh/genomes/renamed_NAU.LB_R.sativus_genome.fasta"
OUTPUT_FASTA="gene_sequences.fa"

# === STEP 1: Extract gene features and convert to BED ===
awk 'BEGIN{FS=OFS="\t"} $3 == "gene" {
    match($9, /ID=([^;]+)/, a)
    id = a[1]
    print $1, $4 - 1, $5, id, ".", $7
}' "$GFF_FILE" > genes.bed

# === STEP 2: Extract sequences with bedtools ===
bedtools getfasta -fi "$FASTA_FILE" -bed genes.bed -s -name > "$OUTPUT_FASTA"

echo "Gene sequences written to $OUTPUT_FASTA"
