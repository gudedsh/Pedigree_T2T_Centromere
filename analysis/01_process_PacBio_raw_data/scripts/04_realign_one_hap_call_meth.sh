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

prefix="${sample}.${hap}.realign"
log="${prefix}.realign_call_meth.log"
unsorted_bam="${prefix}.unsorted.bam"
sorted_bam="${prefix}.sorted.bam"
tempdir="${prefix}.tmp"

mkdir -p "$tempdir"

{
    echo "[$(date)] Re-align phased reads and call CpG methylation"
    echo "Input BAM: $input_bam"
    echo "Sample: $sample"
    echo "Haplotype: $hap"
    echo "Reference: $ref"
    echo "Threads: $threads"

    pbmm2 --version
    samtools --version | head -n 1
    aligned_bam_to_cpg_scores --version

    echo "[$(date)] Running pbmm2 align"

    pbmm2 align \
      "$ref" \
      "$input_bam" \
      "$unsorted_bam" \
      --preset HIFI \
      --sample "${sample}_${hap}" \
      -j "$threads" \
      --log-level INFO

    echo "[$(date)] Sorting BAM"

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

    echo "[$(date)] Calling CpG methylation"

    meth_prefix="${prefix}.MAPQ1.depth1"

    aligned_bam_to_cpg_scores \
      --bam "$sorted_bam" \
      --output-prefix "$meth_prefix" \
      --pileup-mode model \
      --threads "$threads" \
      --min-coverage 1 \
      --min-mapq 1

    echo "[$(date)] Sorting and indexing methylation BED outputs"

    for f in ${meth_prefix}*.bed.gz; do
        [[ -e "$f" ]] || continue

        base=$(basename "$f" .gz)
        sorted_bed="${base}.sorted.bed"

        zcat "$f" | sortBed -i - > "$sorted_bed"
        bgzip -f "$sorted_bed"
        tabix -f -p bed "${sorted_bed}.gz"

        mv -f "${sorted_bed}.gz" "$f"
        mv -f "${sorted_bed}.gz.tbi" "${f}.tbi"
    done

    echo "[$(date)] Done."
    echo "Output BAM: $outdir/$sorted_bam"
    echo "Methylation prefix: $outdir/$meth_prefix"

} 2>&1 | tee "$log"
