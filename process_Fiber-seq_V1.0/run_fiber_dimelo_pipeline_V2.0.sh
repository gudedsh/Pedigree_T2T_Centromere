#!/bin/bash
#SBATCH --job-name=fiber_dimelo
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --output=fiber_dimelo.%j.out
#SBATCH --error=fiber_dimelo.%j.err

set -euo pipefail

########################################
## Usage
########################################

if [[ $# -lt 1 ]]; then
    echo "Usage:"
    echo "  sbatch run_fiber_dimelo_pipeline.sbatch config.sh"
    exit 1
fi

CONFIG=$(realpath "$1")

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: config file not found: $CONFIG"
    exit 1
fi

source "$CONFIG"

########################################
## Basic checks
########################################

: "${SAMPLE:?ERROR: SAMPLE not set}"
: "${INPUT_BAM:?ERROR: INPUT_BAM not set}"
: "${MAT_REF:?ERROR: MAT_REF not set}"
: "${PAT_REF:?ERROR: PAT_REF not set}"
: "${OUTDIR:?ERROR: OUTDIR not set}"

INPUT_BAM=$(realpath "$INPUT_BAM")
MAT_REF=$(realpath "$MAT_REF")
PAT_REF=$(realpath "$PAT_REF")
OUTDIR=$(realpath "$OUTDIR")

MIN_LEN="${MIN_LEN:-1000}"
MG_CUTOFF="${MG_CUTOFF:-0}"

TOTAL_THREADS="${SLURM_CPUS_PER_TASK:-${TOTAL_THREADS:-8}}"
PARALLEL_INIT_JOBS="${PARALLEL_INIT_JOBS:-3}"
PARALLEL_STEP4_JOBS="${PARALLEL_STEP4_JOBS:-2}"

INIT_THREADS=$(( TOTAL_THREADS / PARALLEL_INIT_JOBS ))
STEP4_THREADS=$(( TOTAL_THREADS / PARALLEL_STEP4_JOBS ))

if [[ "$INIT_THREADS" -lt 1 ]]; then INIT_THREADS=1; fi
if [[ "$STEP4_THREADS" -lt 1 ]]; then STEP4_THREADS=1; fi

mkdir -p "$OUTDIR"

LOGDIR="$OUTDIR/logs"
mkdir -p "$LOGDIR"

MAIN_LOG="$LOGDIR/${SAMPLE}.pipeline.log"
VERSION_LOG="$LOGDIR/${SAMPLE}.software_versions.txt"

script_dir=$(dirname "$(realpath "$0")")

{
    echo "========================================"
    echo "Fiber-seq / DiMeLo-seq pipeline"
    echo "Date: $(date)"
    echo "Sample: $SAMPLE"
    echo "Input BAM: $INPUT_BAM"
    echo "Mat ref: $MAT_REF"
    echo "Pat ref: $PAT_REF"
    echo "Outdir: $OUTDIR"
    echo "Total threads: $TOTAL_THREADS"
    echo "Initial parallel jobs: $PARALLEL_INIT_JOBS"
    echo "Threads per initial job: $INIT_THREADS"
    echo "Step4 parallel jobs: $PARALLEL_STEP4_JOBS"
    echo "Threads per Step4 job: $STEP4_THREADS"
    echo "SLURM job ID: ${SLURM_JOB_ID:-NA}"
    echo "SLURM mem: ${SLURM_MEM_PER_NODE:-NA}"
    echo "========================================"
} | tee "$MAIN_LOG"

########################################
## Software check and version record
########################################

check_cmd() {
    cmd="$1"
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: required command not found: $cmd" | tee -a "$MAIN_LOG"
        exit 1
    fi
}

record_version() {
    name="$1"
    shift
    echo "## $name" >> "$VERSION_LOG"
    { "$@" 2>&1 || true; } | head -n 5 >> "$VERSION_LOG"
    echo "" >> "$VERSION_LOG"
}

echo "[$(date)] Checking software..." | tee -a "$MAIN_LOG"

required_cmds=(
    ft
    samtools
    pbmm2
    aligned_bam_to_cpg_scores
    sortBed
    bgzip
    tabix
    awk
    join
    sort
    zcat
)

for cmd in "${required_cmds[@]}"; do
    check_cmd "$cmd"
done

echo "Software version record" > "$VERSION_LOG"
echo "Date: $(date)" >> "$VERSION_LOG"
echo "" >> "$VERSION_LOG"

record_version "ft" ft --version
record_version "samtools" samtools --version
record_version "pbmm2" pbmm2 --version
record_version "aligned_bam_to_cpg_scores" aligned_bam_to_cpg_scores --version
record_version "sortBed" sortBed -version
record_version "bgzip" bgzip --version
record_version "tabix" tabix --version

echo "[$(date)] Software check passed." | tee -a "$MAIN_LOG"

########################################
## Output paths
########################################

STEP1_DIR="$OUTDIR/01_predict_m6a_qc"
STEP2_DIR="$OUTDIR/02_align_unphased"
STEP3_DIR="$OUTDIR/03_phase_reads"
STEP4_DIR="$OUTDIR/04_realign_phased_call_meth"

mkdir -p "$STEP1_DIR" "$STEP2_DIR" "$STEP3_DIR" "$STEP4_DIR"

PRED_BAM="$STEP1_DIR/${SAMPLE}.predict_m6a.bam"

MAT_UNPHASED_BAM="$STEP2_DIR/Mat/${SAMPLE}.Mat.unphased.sorted.bam"
PAT_UNPHASED_BAM="$STEP2_DIR/Pat/${SAMPLE}.Pat.unphased.sorted.bam"

MAT_READS_BAM="$STEP3_DIR/${SAMPLE}.MatReads.bam"
PAT_READS_BAM="$STEP3_DIR/${SAMPLE}.PatReads.bam"

########################################
## Step 1 + Step 2 Mat/Pat in parallel
########################################

echo "[$(date)] Starting Step1 and Step2 Mat/Pat in parallel..." | tee -a "$MAIN_LOG"

bash "$script_dir/scripts/01_predict_m6a_qc.sh" \
    "$INPUT_BAM" \
    "$INIT_THREADS" \
    "$SAMPLE" \
    "$STEP1_DIR" \
    > "$LOGDIR/${SAMPLE}.step1.predict_m6a_qc.log" 2>&1 &
pid_step1=$!

bash "$script_dir/scripts/02_align_one_hap.sh" \
    "$INPUT_BAM" \
    "$SAMPLE" \
    "Mat" \
    "$MAT_REF" \
    "$STEP2_DIR/Mat" \
    "$INIT_THREADS" \
    > "$LOGDIR/${SAMPLE}.step2.Mat.align.log" 2>&1 &
pid_step2_mat=$!

bash "$script_dir/scripts/02_align_one_hap.sh" \
    "$INPUT_BAM" \
    "$SAMPLE" \
    "Pat" \
    "$PAT_REF" \
    "$STEP2_DIR/Pat" \
    "$INIT_THREADS" \
    > "$LOGDIR/${SAMPLE}.step2.Pat.align.log" 2>&1 &
pid_step2_pat=$!

wait "$pid_step1"
wait "$pid_step2_mat"
wait "$pid_step2_pat"

echo "[$(date)] Step1 and Step2 finished." | tee -a "$MAIN_LOG"

########################################
## Check Step1/2 outputs
########################################

for f in "$PRED_BAM" "$MAT_UNPHASED_BAM" "$PAT_UNPHASED_BAM"; do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: expected output missing or empty: $f" | tee -a "$MAIN_LOG"
        exit 1
    fi
done

########################################
## Step 3: phase reads
########################################

echo "[$(date)] Starting Step3 phase reads..." | tee -a "$MAIN_LOG"

bash "$script_dir/scripts/03_phase_reads.sh" \
    "$MAT_UNPHASED_BAM" \
    "$PAT_UNPHASED_BAM" \
    "$PRED_BAM" \
    "$SAMPLE" \
    "$STEP3_DIR" \
    "$TOTAL_THREADS" \
    "$MIN_LEN" \
    "$MG_CUTOFF" \
    > "$LOGDIR/${SAMPLE}.step3.phase_reads.log" 2>&1

echo "[$(date)] Step3 finished." | tee -a "$MAIN_LOG"

for f in "$MAT_READS_BAM" "$PAT_READS_BAM"; do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: expected phased BAM missing or empty: $f" | tee -a "$MAIN_LOG"
        exit 1
    fi
done

########################################
## Step 4 Mat/Pat in parallel
########################################

echo "[$(date)] Starting Step4 Mat/Pat in parallel..." | tee -a "$MAIN_LOG"

bash "$script_dir/scripts/04_realign_one_hap_call_meth.sh" \
    "$MAT_READS_BAM" \
    "$SAMPLE" \
    "Mat" \
    "$MAT_REF" \
    "$STEP4_DIR/Mat" \
    "$STEP4_THREADS" \
    > "$LOGDIR/${SAMPLE}.step4.Mat.realign_call_meth.log" 2>&1 &
pid_step4_mat=$!

bash "$script_dir/scripts/04_realign_one_hap_call_meth.sh" \
    "$PAT_READS_BAM" \
    "$SAMPLE" \
    "Pat" \
    "$PAT_REF" \
    "$STEP4_DIR/Pat" \
    "$STEP4_THREADS" \
    > "$LOGDIR/${SAMPLE}.step4.Pat.realign_call_meth.log" 2>&1 &
pid_step4_pat=$!

wait "$pid_step4_mat"
wait "$pid_step4_pat"

echo "[$(date)] Step4 finished." | tee -a "$MAIN_LOG"

########################################
## Final summary
########################################

echo "========================================" | tee -a "$MAIN_LOG"
echo "Pipeline finished successfully." | tee -a "$MAIN_LOG"
echo "Output directory:" | tee -a "$MAIN_LOG"
echo "  $OUTDIR" | tee -a "$MAIN_LOG"
echo "Software versions:" | tee -a "$MAIN_LOG"
echo "  $VERSION_LOG" | tee -a "$MAIN_LOG"
echo "Main log:" | tee -a "$MAIN_LOG"
echo "  $MAIN_LOG" | tee -a "$MAIN_LOG"
echo "========================================" | tee -a "$MAIN_LOG"
