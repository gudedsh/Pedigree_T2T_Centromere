#!/bin/bash
# Active strict error-trapping switches
set -euo pipefail

if [[ $# -lt 4 ]]; then
    echo "Usage: $0 <input_bam> <threads> <max_reads> <output_bam>"
    exit 1
fi

BAM="$1"
THREADS="$2"
MAX_READS="$3"
OUTPUT_BAM="$4"

PREFIX="${OUTPUT_BAM%.bam}"
HEADER_SAM="${PREFIX}_header.tmp.sam"

echo "[TEST MODE] Subsampling process initiated."
echo "[TEST MODE] Extracting ~1% reads (Hard-capped at a maximum of ${MAX_READS} alignment lines)..."

# Safely isolate the complete BAM header block
samtools view -H "$BAM" > "$HEADER_SAM"

# -s 0.01 selects roughly 1% of reads uniformly across the genome
samtools view -@ "$THREADS" -s 0.01 "$BAM" | head -n "$MAX_READS" >> "$HEADER_SAM"

# Recompile the streaming SAM back into a structured, indexed binary BAM
samtools view -@ "$THREADS" -b "$HEADER_SAM" -o "$OUTPUT_BAM"

# Clean up intermediate header text artifacts
rm -f "$HEADER_SAM"

echo "[TEST MODE] Subsampling successfully finished. Target test file generated: $OUTPUT_BAM"
