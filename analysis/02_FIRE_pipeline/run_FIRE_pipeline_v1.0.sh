#!/bin/bash

# Activate strict error trapping switches
set -euo pipefail

# Print clean usage helper instructions
usage() {
    echo "Usage: bash $0 --sample <prefix> --bam <input.bam> --ref <reference.fa> --outdir <output_dir> [options]"
    echo ""
    echo "Required Arguments:"
    echo "  --sample    Name/Prefix identifier of the genomic lineage track (e.g., PAN010_iPSC_Mat)."
    echo "  --bam       Absolute or relative file path to input coordinates alignment file (.bam)."
    echo "  --ref       Absolute or relative file path to the native genomic assembly (.fa / .fasta)."
    echo "  --outdir    Target workspace output directory."
    echo ""
    echo "Optional Options:"
    echo "  --threads   Number of compute cores allocated for local pipeline steps [Default: 10]."
    echo "  --chrsize   Target custom chromosome text mapping file (.txt) [Default: Auto-derived from reference folder]."
    exit 1
}

# Assign flexible fallback options variables
SAMPLE=""
BAM=""
REF=""
OUTDIR=""
THREADS="10"
CHRSIZE=""

# Parse positional options parameters
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample) SAMPLE="$2"; shift 2 ;;
        --bam)    BAM="$2"; shift 2 ;;
        --ref)    REF="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --chrsize) CHRSIZE="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "[ERROR] Unknown option parameter argument: $1"; usage ;;
    esac
done

# Perform parameter requirement validations
if [ -z "$SAMPLE" ] \vert{}\vert{} [ -z "$BAM" ] || [ -z "$REF" ] \vert{}\vert{} [ -z "$OUTDIR" ]; then
    echo "[ERROR] Missing mandatory structural inputs."
    usage
fi

# Resolve execution variables to solid absolute paths
BAM_ABS=$(realpath "$BAM")
REF_ABS=$(realpath "$REF")
mkdir -p "$OUTDIR"
OUTDIR_ABS=$(realpath "$OUTDIR")

# Core Environment Variable Configurations
export SNAKEMAKE_CONDA_PREFIX=/scratch/devtools/dongs/miniconda3/envs/fibertools
export APPTAINER_CACHEDIR=/wanglab/sdong/miniconda3/envs/snakemake/apptainer-cache
FIRE_BIN="/scratch/devtools/dongs/FIRE/fire"
BIGWIG_TO_BG="/scratch/devtools/dongs/ucsc/bigWigToBedGraph"
BG_TO_BIGWIG="/scratch/devtools/dongs/ucsc/bedGraphToBigWig"

# Auto-derive underlying pedigree metadata signatures
IND=$(echo "$SAMPLE" | grep -o 'PAN[0-9]*' || echo "GENOME")
HAP=$(echo "$SAMPLE" | grep -o '[MP]at' || echo "Hap")

# Locate or set target chromosome sizing metrics files
if [ -z "$CHRSIZE" ]; then
    REF_DIR=$(dirname "$REF_ABS")
    CHRSIZE_MATCH=$(ls${REF_DIR}/*${IND}*${HAP}*_chrsize.txt 2>/dev/null | head -n 1 || true)
    if [ -f "$CHRSIZE_MATCH" ]; then
        CHRSIZE_ABS="$CHRSIZE_MATCH"
    else
        # Fallback tracking if specific multi-string file names do not match template
        CHRSIZE_ABS="${REF_DIR}/${IND}_assembly.v1.0_${HAP}_chr1-22-XY_chrsize.txt"
    fi
else
    CHRSIZE_ABS=$(realpath "$CHRSIZE")
fi

# Initialize logging channels
LOG_FILE="${OUTDIR_ABS}/${SAMPLE}_pipeline_execution.log"
echo "=========================================================" | tee -a "$LOG_FILE"
echo " [INFO] Initializing FIRE Pipeline Wrapper Step" | tee -a "$LOG_FILE"
echo " [TIME] $(date)" \vert{} tee -a "$LOG_FILE"
echo "=========================================================" | tee -a "$LOG_FILE"
echo ">> Sample Prefix : $SAMPLE" \vert{} tee -a "$LOG_FILE"
echo ">> Input BAM     : $BAM_ABS" \vert{} tee -a "$LOG_FILE"
echo ">> Reference FA  : $REF_ABS" \vert{} tee -a "$LOG_FILE"
echo ">> Size Blueprint: $CHRSIZE_ABS" \vert{} tee -a "$LOG_FILE"
echo ">> Threads Count : $THREADS" \vert{} tee -a "$LOG_FILE"
echo ">> Output Workspace: $OUTDIR_ABS" \vert{} tee -a "$LOG_FILE"

# 1. Coordinate and verify indexing benchmarks
echo "[INFO] Running samtools validation benchmarks..." | tee -a "$LOG_FILE"
if [ ! -f "${REF_ABS}.fai" ]; then
    echo "[WARN] Reference index missing. Generating .fai..." | tee -a "$LOG_FILE"
    samtools faidx "$REF_ABS" -@ "$THREADS"
fi
if [ ! -f "${BAM_ABS}.bai" ]; then
    echo "[WARN] BAM index missing. Generating .bai..." | tee -a "$LOG_FILE"
    samtools index -b "$BAM_ABS" -@ "$THREADS"
fi

# 2. Build tracking configuration workspaces
CONFIG_DIR="${OUTDIR_ABS}/config"
mkdir -p "$CONFIG_DIR"
MANIFEST_TBL="${CONFIG_DIR}/${SAMPLE}.tbl"
CONFIG_YAML="${CONFIG_DIR}/${SAMPLE}.yaml"

echo -e "sample\tbam\n${SAMPLE}\t${BAM_ABS}" > "$MANIFEST_TBL"
echo -e "ref: ${REF_ABS}\nref_name: ${IND}_${HAP}\nmanifest: ${MANIFEST_TBL}" > "$CONFIG_YAML"

# 3. Deploy the core FIRE calculation pipeline
echo "[INFO] Launching underlying core FIRE computational framework..." | tee -a "$LOG_FILE"
cd "$OUTDIR_ABS"

if ${FIRE_BIN} --configfile "$CONFIG_YAML" >> "$LOG_FILE" 2>&1; then
    echo "[SUCCESS] Core FIRE computational process finished successfully." | tee -a "$LOG_FILE"
else
    echo "[ERROR] FIRE run crashed. Review full execution metrics inside: $LOG_FILE"
    exit 1
fi


echo "=========================================================" | tee -a "$LOG_FILE"
echo " [SUCCESS] Entire Module Wrapper Process Fully Completed." | tee -a "$LOG_FILE"
echo " [TIME] $(date)" \vert{} tee -a "$LOG_FILE"
echo "=========================================================" | tee -a "$LOG_FILE"
