#!/bin/bash

# Usage: ./generate_hg38_5kb_bins.sh /path/to/hg38.fa

set -e

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 /path/to/hg38.fa"
    exit 1
fi

FASTA="$1"
BASENAME=$(basename "$FASTA" .fa)
GENOME="${BASENAME}.genome"
GENOME_FILTERED="${BASENAME}.genome.filtered"
BED="${BASENAME}_5kb_bins.bed"

# Step 1: Index genome fasta using samtools env
source /home/mattia/miniconda3/etc/profile.d/conda.sh
conda activate samtools
samtools faidx "$FASTA"

# Step 2: Create genome file with chrom sizes
cut -f1,2 "${FASTA}.fai" > "$GENOME"

# Step 3: Filter for canonical chromosomes (chr1–chr22, chrX, chrY)
egrep -w "chr([1-9]|1[0-9]|2[0-2]|X|Y|M)" "$GENOME" > "$GENOME_FILTERED"

# Step 4: Generate 5kb bins as BED file using bedtools env
conda activate bedtools
bedtools makewindows -g "$GENOME_FILTERED" -w 5000 > "$BED"

echo "Generated $BED (canonical chromosomes only)"
