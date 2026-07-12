#!/bin/bash

# Activate strict error trapping switches
set -euo pipefail

# Print clean usage helper instructions
usage() {
    echo "Usage: bash $0 --tumor-bam <tumor.bam> --normal-bam <normal.bam> --ref <reference.fa> --outdir <output_dir> [options]"
    echo ""
    echo "Required Arguments:"
    echo "  --tumor-bam   Absolute or relative file path to the Tumor input track (.bam)."
    echo "  --normal-bam  Absolute or relative file path to the Normal/Control matched input track (.bam)."
    echo "  --ref         Absolute or relative file path to the individual assembly reference (.fa / .fasta)."
    echo "  --outdir      Target workspace output directory where VCFs and logs will be saved."
    echo ""
    echo "Optional Options:"
    echo "  --threads     Number of parallel Native PairHMM compute threads allocated [Default: 10]."
    exit 1
}

# Assign flexible fallback options variables
TUMOR_BAM=""
NORMAL_BAM=""
REF=""
OUTDIR=""
THREADS="10"

# Parse positional options parameters
while [[ $# -gt 0 ]]; do
    case "$1" in
        --tumor-bam)  TUMOR_BAM="$2"; shift 2 ;;
        --normal-bam) NORMAL_BAM="$2"; shift 2 ;;
        --ref)        REF="$2"; shift 2 ;;
        --outdir)     OUTDIR="$2"; shift 2 ;;
        --threads)    THREADS="$2"; shift 2 ;;
        -h|--help)    usage ;;
        *) echo "[ERROR] Unknown option parameter argument: $1"; usage ;;
    esac
done

# Perform parameter requirement validations
if [ -z "$TUMOR_BAM" ] \vert{}\vert{} [ -z "$NORMAL_BAM" ] || [ -z "$REF" ] \vert{}\vert{} [ -z "$OUTDIR" ]; then
    echo "[ERROR] Missing mandatory structural input parameters."
    usage
fi

# Resolve execution tracks to concrete absolute paths
TUMOR_ABS=$(realpath "$TUMOR_BAM")
NORMAL_ABS=$(realpath "$NORMAL_BAM")
REF_ABS=$(realpath "$REF")
mkdir -p "$OUTDIR"
OUTDIR_ABS=$(realpath "$OUTDIR")

# Deriving target descriptors for output naming
NAME=$(basename "$TUMOR_ABS" .bam)
FUNC="gatk_Mutect2_local"

# Establish destination file tokens
PREFIX="${OUTDIR_ABS}/${NAME}"
RAW_VCF="${PREFIX}_${FUNC}.vcf.gz"
F1R2_TAR="${PREFIX}_f1r2.tar.gz"
ORIEN_MODEL="${PREFIX}_read-orientation-model.tar.gz"
FILTERED_VCF="${PREFIX}_${FUNC}_filtered.vcf.gz"

echo "========================================================="
echo " [INFO] Initializing GATK4 Mutect2 Somatic Call Pipeline"
echo " [TIME] $(date)"
echo "========================================================="
echo ">> Target Tumor BAM : $TUMOR_ABS"
echo ">> Matched Normal   : $NORMAL_ABS"
echo ">> Reference FASTA  : $REF_ABS"
echo ">> PairHMM Threads  : $THREADS Cores"
echo ">> Output Workspace : $OUTDIR_ABS"

# 1. Verification of sequence indexing benchmarks
echo "[INFO] Running indexing benchmark checks..."
if [ ! -f "${REF_ABS}.fai" ]; then
    echo "[WARN] Reference index missing. Generating .fai..."
    samtools faidx "$REF_ABS" -@ "$THREADS"
fi
if [ ! -f "${TUMOR_ABS}.bai" ]; then
    echo "[WARN] Tumor BAM index missing. Generating .bai..."
    samtools index -b "$TUMOR_ABS" -@ "$THREADS"
fi
if [ ! -f "${NORMAL_ABS}.bai" ]; then
    echo "[WARN] Normal BAM index missing. Generating .bai..."
    samtools index -b "$NORMAL_ABS" -@ "$THREADS"
fi

# 2. Automatically parse Read Group Sample Names (SM)
echo "[INFO] Extracting sample tokens from alignment metadata headers..."
TUMOR_SM=$(samtools view -H "$TUMOR_ABS" | grep '^@RG' | awk '{for(i=1;i<=NF;i++) if($i~/^SM:/) print$i}' | sed 's/SM://g' | sort -u | head -n 1)
NORMAL_SM=$(samtools view -H "$NORMAL_ABS" | grep '^@RG' | awk '{for(i=1;i<=NF;i++) if($i~/^SM:/) print$i}' | sed 's/SM://g' | sort -u | head -n 1)

if [ -z "$TUMOR_SM" ] \vert{}\vert{} [ -z "$NORMAL_SM" ]; then
    echo "[ERROR] Failed to extract valid SM tags from Read Group headers. Ensure BAMs contain proper @RG tags."
    exit 1
fi

echo " -> Detected Tumor SM Tag  : $TUMOR_SM"
echo " -> Detected Normal SM Tag : $NORMAL_SM"

# Change directory into the output folder for GATK logs collection
cd "$OUTDIR_ABS"

# 3. Step 01: Core Mutect2 Calling
echo "[INFO] [STAGE 1/3] Running GATK Mutect2..."
if gatk --java-options "-Xmx32g" Mutect2 \
  -R "$REF_ABS" \
  -I "$TUMOR_ABS" -tumor "$TUMOR_SM" \
  -I "$NORMAL_ABS" -normal "$NORMAL_SM" \
  --min-base-quality-score 20 \
  --minimum-mapping-quality 20 \
  --native-pair-hmm-threads "$THREADS" \
  --max-reads-per-alignment-start 0 \
  -O "$RAW_VCF" \
  --f1r2-tar-gz "$F1R2_TAR" &> "01_${NAME}_gatk_Mutect2.log"; then
    echo " -> Mutect2 core variant calling successfully completed."
else
    echo "[ERROR] Mutect2 execution failed. Check log: 01_${NAME}_gatk_Mutect2.log"
    exit 1
    _
fi

# 4. Step 02: Learn Read Orientation Model
echo "[INFO] [STAGE 2/3] Running GATK LearnReadOrientationModel..."
if gatk LearnReadOrientationModel \
  -I "$F1R2_TAR" \
  -O "$ORIEN_MODEL" &> "02_${NAME}_gatk_LearnReadOrientationModel.log"; then
    echo " -> Bias model generation successfully completed."
else
    echo "[ERROR] LearnReadOrientationModel failed. Check log: 02_${NAME}_gatk_LearnReadOrientationModel.log"
    exit 1
fi

# 5. Step 03: Apply Somatic Filters
echo "[INFO] [STAGE 3/3] Running GATK FilterMutectCalls..."
if gatk FilterMutectCalls \
  -R "$REF_ABS" \
  -V "$RAW_VCF" \
  --ob-priors "$ORIEN_MODEL" \
  -O "$FILTERED_VCF" &> "03_${NAME}_gatk_FilterMutectCalls.log"; then
    echo "========================================================="
    echo " [SUCCESS] GATK Mutect2 Somatic Calling Fully Completed."
    echo " [FINAL OUTPUT] $FILTERED_VCF"
    echo "========================================================="
else
    echo "[ERROR] FilterMutectCalls execution failed. Check log: 03_${NAME}_gatk_FilterMutectCalls.log"
    exit 1
fi
