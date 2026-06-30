#!/bin/bash
#SBATCH --job-name=call_FAS_from_FIRE
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G

set -euo pipefail

usage() {
cat <<EOF
Usage:
  bash call_FAS_pipeline_V1.0.sh \\
    -fa fa_files.tsv \\
    -cen active_hor_cdr.bed \\
    -o call_FAS_output \\
    -bin 2000 \\
    -s 5 \\
    -exp 0 \\
    -out_CDR_cutoff 0.10 \\
    -in_CDR_cutoff 0.05 \\
    -t 8

Required:
  -fa      two-column TSV: sample_id coverage_bed_gz
  -cen     BED-like file: chr activeHOR_start activeHOR_end cdr_start cdr_end

Optional:
  -o       output directory [default: call_FAS_output]
  -bin     bin size [default: 2000]
  -s       smoothing window size [default: 5]
  -exp     symmetric expansion size for CDR only [default: 0]
  -out_CDR_cutoff  winsorization fraction outside CDR [default: 0.10]
  -in_CDR_cutoff   winsorization fraction inside CDR [default: 0.05]
  -t       threads [default: SLURM_CPUS_PER_TASK or 8]
EOF
}

fa_file=""
activehor=""
outdir="call_FAS_output"
winsize=5000
smooth_window=10
expansion_size=0
out_CDR_cutoff=0.10
in_CDR_cutoff=0.05
threads="${SLURM_CPUS_PER_TASK:-8}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -fa) fa_file="$2"; shift 2 ;;
    -cen) activehor="$2"; shift 2 ;;
    -o) outdir="$2"; shift 2 ;;
    -bin) winsize="$2"; shift 2 ;;
    -s) smooth_window="$2"; shift 2 ;;
    -exp) expansion_size="$2"; shift 2 ;;
    -out_CDR_cutoff) out_CDR_cutoff="$2"; shift 2 ;;
    -in_CDR_cutoff) in_CDR_cutoff="$2"; shift 2 ;;
    -t) threads="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

[[ -z "$fa_file" || -z "$activehor" ]] && { usage; exit 1; }

[[ -f "$fa_file" ]] || { echo "ERROR: fa file not found: $fa_file"; exit 1; }
[[ -f "$activehor" ]] || { echo "ERROR: active HOR file not found: $activehor"; exit 1; }

command -v tabix >/dev/null || { echo "ERROR: tabix not found"; exit 1; }
command -v Rscript >/dev/null || { echo "ERROR: Rscript not found"; exit 1; }

outdir=$(realpath -m "$outdir")
mkdir -p "$outdir"

workdir="${outdir}/intermediate_files"
extract_dir="${workdir}/activeHOR_coverage"
plotdir="${outdir}/FAS_QC_plots"

mkdir -p "$workdir" "$extract_dir" "$plotdir"

log="${outdir}/call_FAS.log"

{
  echo "=========== Call FAS from FIRE coverage ============"
  date
  echo "FA file: $fa_file"
  echo "Active HOR / CDR file: $activehor"
  echo "Output dir: $outdir"
  echo "Bin size: $winsize"
  echo "Smooth window: $smooth_window"
  echo "CDR expansion: $expansion_size"
  echo "out_CDR_cutoff: $out_CDR_cutoff"
  echo "in_CDR_cutoff: $in_CDR_cutoff"
  echo "Threads: $threads"
} > "$log"

active_clean="${workdir}/active_hor_cdr.clean.bed"

awk 'BEGIN{OFS="\t"}
  /^#/ || NF==0 {next}
  NF < 5 {
    print "ERROR: active HOR file must have at least 5 columns at line " NR > "/dev/stderr"
    exit 1
  }
  {
    print $1,$2,$3,$4,$5
  }
' "$activehor" > "$active_clean"

combined_active_cov="${outdir}/FAS_activeHOR_coverage.tsv"
echo -e "sample\tchr\tstart\tend\tfire\tlinker\tnucleosome" > "$combined_active_cov"

echo "[$(date)] Extracting active HOR coverage by tabix..." | tee -a "$log"

while read -r sample cov; do
  [[ -z "${sample:-}" ]] && continue
  [[ "$sample" =~ ^# ]] && continue

  [[ -f "$cov" ]] || { echo "ERROR: coverage file not found: $cov"; exit 1; }

  if [[ ! -f "${cov}.tbi" ]]; then
    echo "ERROR: tabix index not found for $cov" | tee -a "$log"
    echo "Please prepare with:" | tee -a "$log"
    echo "  sort -k1,1 -k2,2n input.bed | bgzip > input.sorted.bed.gz" | tee -a "$log"
    echo "  tabix -p bed input.sorted.bed.gz" | tee -a "$log"
    exit 1
  fi

  sample_out="${extract_dir}/${sample}.activeHOR.coverage.tsv"
  : > "$sample_out"

  echo "  sample: $sample" | tee -a "$log"

  while read -r chr active_start active_end cdr_start cdr_end; do
    tabix "$cov" "${chr}:${active_start}-${active_end}" 2>/dev/null || true
  done < "$active_clean" | \
  awk -v s="$sample" 'BEGIN{OFS="\t"} NF>=6 {print s,$1,$2,$3,$4,$5,$6}' \
  > "$sample_out"

  cat "$sample_out" >> "$combined_active_cov"

done < "$fa_file"

echo "[$(date)] Active HOR coverage written to:" | tee -a "$log"
echo "  $combined_active_cov" | tee -a "$log"

rcode=$(realpath call_fas_from_fire_coverage.r)

[[ -f "$rcode" ]] || {
  echo "ERROR: R code not found: $rcode" | tee -a "$log"
  exit 1
}

fas_out="${outdir}/FAS_results.tsv"
bin_out="${outdir}/FAS_bin_summary.tsv"

echo "[$(date)] Calling FAS in R..." | tee -a "$log"

Rscript "$rcode" \
  --coverage "$combined_active_cov" \
  --cen "$active_clean" \
  --out "$fas_out" \
  --bin_out "$bin_out" \
  --plotdir "$plotdir" \
  --bin_size "$winsize" \
  --smooth_k "$smooth_window" \
  --cdr_expansion "$expansion_size" \
  --out_CDR_cutoff "$out_CDR_cutoff" \
  --in_CDR_cutoff "$in_CDR_cutoff" \
  >> "$log" 2>&1

echo "[$(date)] Done." | tee -a "$log"
echo "FAS result: $fas_out"
echo "Bin summary: $bin_out"
echo "QC plots: $plotdir"
