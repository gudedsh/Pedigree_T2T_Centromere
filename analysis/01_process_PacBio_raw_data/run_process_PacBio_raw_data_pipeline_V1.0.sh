#!/bin/bash
set -euo pipefail

usage() {
cat <<EOF
Usage:
  bash run_fiber_dimelo_pipeline.sh \\
    --input-bam input.bam \\
    --sample sample_name \\
    --mat-ref Mat.fa \\
    --pat-ref Pat.fa \\
    --outdir output_dir \\
    --threads 8

Required:
  --input-bam   raw PacBio BAM
  --sample      sample name
  --mat-ref     maternal haplotype fasta
  --pat-ref     paternal haplotype fasta
  --outdir      output directory

Optional:
  --threads     threads [default: 8]
  --min-len     minimum read length for phasing [default: 1000]
  --mg-cutoff   mg difference cutoff [default: 0]
EOF
}

threads=8
min_len=1000
mg_cutoff=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-bam) input_bam=$(realpath "$2"); shift 2 ;;
    --sample) sample="$2"; shift 2 ;;
    --mat-ref) mat_ref=$(realpath "$2"); shift 2 ;;
    --pat-ref) pat_ref=$(realpath "$2"); shift 2 ;;
    --outdir) outdir=$(realpath "$2"); shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --min-len) min_len="$2"; shift 2 ;;
    --mg-cutoff) mg_cutoff="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

: "${input_bam:?ERROR: --input-bam required}"
: "${sample:?ERROR: --sample required}"
: "${mat_ref:?ERROR: --mat-ref required}"
: "${pat_ref:?ERROR: --pat-ref required}"
: "${outdir:?ERROR: --outdir required}"

script_dir=$(dirname "$(realpath "$0")")
mkdir -p "$outdir"

echo "========================================"
echo "Fiber-seq / DiMeLo-seq integrated pipeline"
echo "Sample:   $sample"
echo "Input:    $input_bam"
echo "Mat ref:  $mat_ref"
echo "Pat ref:  $pat_ref"
echo "Outdir:   $outdir"
echo "Threads:  $threads"
echo "========================================"

########################################
## Step 1. predict m6A and QC
########################################

mkdir -p "$outdir/01_predict_m6a_qc"

bash "$script_dir/scripts/01_predict_m6a_qc.sh" \
  "$input_bam" \
  "$threads" \
  "$sample" \
  "$outdir/01_predict_m6a_qc"

pred_bam="$outdir/01_predict_m6a_qc/${sample}.predict_m6a.bam"

########################################
## Step 2. align unphased reads to Mat and Pat
########################################

mkdir -p "$outdir/02_align_unphased"

bash "$script_dir/scripts/02_align_unphased_to_haps.sh" \
  "$pred_bam" \
  "$sample" \
  "$mat_ref" \
  "$pat_ref" \
  "$outdir/02_align_unphased" \
  "$threads"

mat_bam="$outdir/02_align_unphased/Mat/${sample}.Mat.unphased.sorted.bam"
pat_bam="$outdir/02_align_unphased/Pat/${sample}.Pat.unphased.sorted.bam"

########################################
## Step 3. extract scores and phase reads
########################################

mkdir -p "$outdir/03_phase_reads"

bash "$script_dir/scripts/03_phase_reads.sh" \
  "$mat_bam" \
  "$pat_bam" \
  "$pred_bam" \
  "$sample" \
  "$outdir/03_phase_reads" \
  "$threads" \
  "$min_len" \
  "$mg_cutoff"

mat_reads_bam="$outdir/03_phase_reads/${sample}.MatReads.bam"
pat_reads_bam="$outdir/03_phase_reads/${sample}.PatReads.bam"

########################################
## Step 4. realign phased reads and call methylation
########################################

mkdir -p "$outdir/04_realign_phased_call_meth"

bash "$script_dir/scripts/04_realign_phased_call_meth.sh" \
  "$mat_reads_bam" \
  "$mat_ref" \
  "${sample}_Mat" \
  "$outdir/04_realign_phased_call_meth/Mat" \
  "$threads"

bash "$script_dir/scripts/04_realign_phased_call_meth.sh" \
  "$pat_reads_bam" \
  "$pat_ref" \
  "${sample}_Pat" \
  "$outdir/04_realign_phased_call_meth/Pat" \
  "$threads"

echo "========================================"
echo "Pipeline finished successfully."
echo "Final output:"
echo "  $outdir"
echo "========================================"
