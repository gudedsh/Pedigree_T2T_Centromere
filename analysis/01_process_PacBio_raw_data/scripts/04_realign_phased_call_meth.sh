#!/bin/bash
set -euo pipefail

input_bam=$(realpath "$1")
ref=$(realpath "$2")
sample="$3"
outdir=$(realpath "$4")
threads="${5:-8}"

mkdir -p "$outdir"
cd "$outdir"

prefix="${sample}.realign"
log="${prefix}.log"

{
    echo "[$(date)] Re-align phased reads and call methylation"
    echo "Input BAM: $input_bam"
    echo "Reference: $ref"
    echo "Sample: $sample"
    echo "Threads: $threads"

    pbmm2 --version

    tempdir="${prefix}.tmp"
    mkdir -p "$tempdir"

    pbmm2 align \
      "$ref" \
      "$input_bam" \
      "${prefix}.unsorted.bam" \
      --preset HIFI \
      --sample "$sample" \
      -j "$threads" \
      --log-level INFO

    samtools sort \
      -@ "$threads" \
      -T "$tempdir" \
      -o "${prefix}.sorted.bam" \
      "${prefix}.unsorted.bam"

    samtools index -@ "$threads" "${prefix}.sorted.bam"

    samtools stats -@ "$threads" "${prefix}.sorted.bam" > "${prefix}.sorted.bam.stats.txt"

    rm -rf "$tempdir" "${prefix}.unsorted.bam"

    echo "[$(date)] Calling CpG methylation"

    aligned_bam_to_cpg_scores --version

    aligned_bam_to_cpg_scores \
      --bam "${prefix}.sorted.bam" \
      --output-prefix "${prefix}.MAPQ1.depth1" \
      --pileup-mode model \
      --threads "$threads" \
      --min-coverage 1 \
      --min-mapq 1

    for f in ${prefix}.MAPQ1.depth1*.bed.gz; do
        [[ -e "$f" ]] || continue

        base=$(basename "$f" .gz)
        sorted_bed="${base}.sorted.bed"

        zcat "$f" | sortBed -i - > "$sorted_bed"
        bgzip -f "$sorted_bed"
        tabix -f -p bed "${sorted_bed}.gz"

        mv "${sorted_bed}.gz" "$f"
        mv "${sorted_bed}.gz.tbi" "${f}.tbi"
    done

    echo "[$(date)] Step 4 done."

} 2>&1 | tee "$log"
