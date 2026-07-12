#!/bin/bash

# Activate strict error trapping switches
set -euo pipefail

# Print clean usage helper instructions
usage() {
    echo "Usage: bash $0 --bam <input.bam> --ref <reference.fa> --outdir <output_dir> [options]"
    echo ""
    echo "Required Arguments:"
    echo "  --bam       Absolute or relative file path to the input alignment track (.bam)."
    echo "  --ref       Absolute or relative file path to the individual assembly reference (.fa / .fasta)."
    echo "  --outdir    Target workspace output directory where VCFs and logs will be saved."
    echo ""
    echo "Optional Options:"
    echo "  --threads   Number of parallel shards allocated for DeepVariant [Default: 10]."
    exit 1
}

# Assign flexible fallback options variables
BAM=""
REF=""
OUTDIR=""
THREADS="10"

# Parse positional options parameters
while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam)     BAM="$2"; shift 2 ;;
        --ref)     REF="$2"; shift 2 ;;
        --outdir)  OUTDIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "[ERROR] Unknown option parameter argument: $1"; usage ;;
    esac
done

# Perform parameter requirement validations
if [ -z "$BAM" ] || [ -z "$REF" ] \vert{}\vert{} [ -z "$OUTDIR" ]; then
    echo "[ERROR] Missing mandatory structural input parameters."
    usage
fi

# Resolve execution tracks to concrete absolute paths
BAM_ABS=$(realpath "$BAM")
REF_ABS=$(realpath "$REF")
mkdir -p "$OUTDIR"
OUTDIR_ABS=$(realpath "$OUTDIR")

# Deriving target descriptors for clean downstream logging nomenclature
NAME=$(basename "$BAM_ABS" .bam)
FUNC="deepvariants_callGV"

# Establish target output track paths
OUTPUT_VCF="${OUTDIR_ABS}/${NAME}_${FUNC}.vcf.gz"
OUTPUT_GVCF="${OUTDIR_ABS}/${NAME}_${FUNC}.gvcf.gz"
LOG_FILE="${OUTDIR_ABS}/${NAME}_${FUNC}.run.log"

echo "========================================================="
echo " [INFO] Initializing DeepVariant Variant Calling Workflow"
echo " [TIME] $(date)"
echo "========================================================="
echo ">> Target Input BAM: $BAM_ABS"
echo ">> Reference FASTA : $REF_ABS"
echo ">> Allocated Shards: $THREADS Cores"
echo ">> Output Workspace: $OUTDIR_ABS"
echo ">> Pipeline Log    : $LOG_FILE"

# 1. Verification of sequence indexing benchmarks
echo "[INFO] Running indexing benchmark checks..."
if [ ! -f "${REF_ABS}.fai" ]; then
    echo "[WARN] Reference index missing. Generating .fai natively..."
    samtools faidx "$REF_ABS" -@ "$THREADS"
fi
if [ ! -f "${BAM_ABS}.bai" ]; then
    echo "[WARN] BAM alignment index missing. Generating .bai alignment index..."
    samtools index -b "$BAM_ABS" -@ "$THREADS"
fi

# 2. Deploy native DeepVariant framework calculation engine
echo "[INFO] Launching run_deepvariant internal engine..."
echo "-> Processing details are being streamed to: $LOG_FILE"

if run_deepvariant \
  --model_type=PACBIO \
  --ref="${REF_ABS}" \
  --reads="${BAM_ABS}" \
  --output_vcf="${OUTPUT_VCF}" \
  --output_gvcf="${OUTPUT_GVCF}" \
  --num_shards="${THREADS}" >> "$LOG_FILE" 2>&1; then
    
    echo "========================================================="
    echo " [SUCCESS] DeepVariant Local Process Completed Successfully."
    echo " [OUTPUT 1] VCF:  $OUTPUT_VCF"
    echo " [OUTPUT 2] gVCF: $OUTPUT_GVCF"
    echo "========================================================="
else
    echo "[ERROR] DeepVariant process failed. Review full logs inside: $LOG_FILE"
    exit 1
fi
