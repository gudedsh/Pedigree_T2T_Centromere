#!/bin/bash
set -euo pipefail

input_bam=$(realpath "$1")
sample="$2"
hap="$3"          # Mat or Pat
ref=$(realpath "$4")
outdir=$(realpath "$5")
threads="${6:-8}"

mkdir -p "$outdir"
cd "$outdir"

prefix="${sample}.${hap}.unphased"
log="${prefix}.pbmm2.log"
unsorted_bam="${prefix}.unsorted.bam"
sorted_bam="${prefix}.sorted.bam"
tempdir="${prefix}.tmp"

mkdir -p "$tempdir"

{
    echo "[$(date)] Aligning unphased reads"
    echo "Input BAM: $input_bam"
    echo "Sample: $sample"
    echo "Haplotype: $hap"
    echo "Reference: $ref"
    echo "Threads: $threads"

    pbmm2 --version
    samtools --version | head -n 1

    pbmm2 align \
      "$ref" \
      "$input_bam" \
      "$unsorted_bam" \
      --preset HIFI \
      --sample "$sample" \
      -j "$threads" \
      --log-level INFO

    samtools sort \
      -@ "$threads" \
      -T "$tempdir" \
      -o "$sorted_bam" \
      "$unsorted_bam"

    samtools index -@ "$threads" "$sorted_bam"

    samtools stats \
      -@ "$threads" \
      "$sorted_bam" \
      > "${sorted_bam}.stats.txt"

    rm -rf "$tempdir" "$unsorted_bam"

    echo "[$(date)] Done."
    echo "Output BAM: $outdir/$sorted_bam"

} 2>&1 | tee "$log"
