#!/bin/bash

input_bam=$(realpath "$1")
sample="$2"
mat_ref=$(realpath "$3")
pat_ref=$(realpath "$4")
outdir=$(realpath "$5")
threads="${6:-8}"

mkdir -p "$outdir"

align_one_hap() {
    hap="$1"
    ref="$2"

    hap_dir="$outdir/$hap"
    mkdir -p "$hap_dir"
    cd "$hap_dir"

    prefix="${sample}.${hap}.unphased"
    log="${prefix}.pbmm2.log"
    unsorted_bam="${prefix}.unsorted.bam"
    sorted_bam="${prefix}.sorted.bam"
    tempdir="${prefix}.tmp"

    mkdir -p "$tempdir"

    {
        echo "[$(date)] Aligning unphased reads to $hap reference"
        echo "Input BAM: $input_bam"
        echo "Reference: $ref"
        echo "Threads: $threads"

        pbmm2 --version

        pbmm2 align \
          "$ref" \
          "$input_bam" \
          "$unsorted_bam" \
          --preset HIFI \
          --sample "$sample" \
          -j "$threads" \
          --log-level INFO

        samtools sort -@ "$threads" -T "$tempdir" -o "$sorted_bam" "$unsorted_bam"
        samtools index -@ "$threads" "$sorted_bam"
        samtools stats -@ "$threads" "$sorted_bam" > "${sorted_bam}.stats.txt"

        rm -rf "$tempdir" "$unsorted_bam"

        echo "[$(date)] Done $hap alignment"
    } 2>&1 | tee "$log"
}

align_one_hap "Mat" "$mat_ref"
align_one_hap "Pat" "$pat_ref"

