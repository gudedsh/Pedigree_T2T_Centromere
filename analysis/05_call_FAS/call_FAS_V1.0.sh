#!/usr/bin/env bash
#SBATCH --job-name=call_FAS_from_FIRE
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G

set -Eeuo pipefail

readonly PROGRAM_NAME="$(basename "$0")"
readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"

usage() {
  cat <<'EOF'
Usage:
  call_FAS.sh --fiberseq-list FILE --centromere FILE [options]

Required:
  --fiberseq-list, -fa FILE      Two-column TSV: sample_id coverage_bed_gz
  --centromere, -cen FILE        Five-column active-HOR/CDR annotation

Options (defaults reproduce v1.0):
  --output-dir, -o DIR           Output directory [call_FAS_output]
  --bin-size, -bin INT           Bin size in bp [5000]
  --smooth-window, -s INT        Smoothing window in bins [10]
  --cdr-expansion, -exp INT      Symmetric CDR expansion in bp [0]
  --outside-cdr-cutoff, -out_CDR_cutoff FRACTION
                                    Winsorization fraction outside CDR [0.10]
  --inside-cdr-cutoff, -in_CDR_cutoff FRACTION
                                    Winsorization fraction inside CDR [0.05]
  --threads, -t INT              Reserved for v1.0 command compatibility
  -h, --help                     Show this help message

Example:
  bash call_FAS.sh \
    --fiberseq-list example_data/fs_files.tsv \
    --centromere example_data/activeHOR_CDR.tsv \
    --output-dir example_output
EOF
}

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
require_value() { [[ $# -ge 2 && -n "${2:-}" ]] || die "Option '$1' requires a value."; }
require_command() { command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"; }
is_positive_integer() { [[ "$1" =~ ^[1-9][0-9]*$ ]]; }
is_integer() { [[ "$1" =~ ^-?[0-9]+$ ]]; }
is_fraction() {
  [[ "$1" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]] &&
    awk -v x="$1" 'BEGIN { exit !(x < 0.5) }'
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

while (($#)); do
  case "$1" in
    --fiberseq-list|-fa) require_value "$@"; fa_file="$2"; shift 2 ;;
    --centromere|-cen) require_value "$@"; activehor="$2"; shift 2 ;;
    --output-dir|-o) require_value "$@"; outdir="$2"; shift 2 ;;
    --bin-size|-bin) require_value "$@"; winsize="$2"; shift 2 ;;
    --smooth-window|-s) require_value "$@"; smooth_window="$2"; shift 2 ;;
    --cdr-expansion|-exp) require_value "$@"; expansion_size="$2"; shift 2 ;;
    --outside-cdr-cutoff|-out_CDR_cutoff) require_value "$@"; out_CDR_cutoff="$2"; shift 2 ;;
    --inside-cdr-cutoff|-in_CDR_cutoff) require_value "$@"; in_CDR_cutoff="$2"; shift 2 ;;
    --threads|-t) require_value "$@"; threads="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1 (run '${PROGRAM_NAME} --help')" ;;
  esac
done

[[ -n "$fa_file" ]] || die "--fiberseq-list is required."
[[ -n "$activehor" ]] || die "--centromere is required."
[[ -r "$fa_file" ]] || die "Fiber-seq manifest not found or unreadable: $fa_file"
[[ -r "$activehor" ]] || die "Active-HOR/CDR annotation not found or unreadable: $activehor"

is_positive_integer "$winsize" || die "--bin-size must be a positive integer."
is_positive_integer "$smooth_window" || die "--smooth-window must be a positive integer."
is_integer "$expansion_size" || die "--cdr-expansion must be an integer."
is_fraction "$out_CDR_cutoff" || die "--outside-cdr-cutoff must be in [0, 0.5)."
is_fraction "$in_CDR_cutoff" || die "--inside-cdr-cutoff must be in [0, 0.5)."
is_positive_integer "$threads" || die "--threads must be a positive integer."

for tool in awk cut Rscript sort tabix; do require_command "$tool"; done

r_code="${SCRIPT_DIR}/r_code/call_fas_from_fire_coverage.r"
[[ -r "$r_code" ]] || die "R program not found: $r_code"

mkdir -p -- "$outdir"
outdir="$(cd -- "$outdir" && pwd -P)"
workdir="${outdir}/intermediate_files"
extract_dir="${workdir}/activeHOR_coverage"
plotdir="${outdir}/FAS_QC_plots"
mkdir -p -- "$workdir" "$extract_dir" "$plotdir"
log_file="${outdir}/call_FAS.log"

active_clean="${workdir}/active_HOR_CDR.validated.tsv"
awk '
  BEGIN { OFS="\t" }
  /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
  NF < 5 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $4 !~ /^[0-9]+$/ || $5 !~ /^[0-9]+$/ {
    printf "ERROR: invalid five-column annotation at line %d.\n", NR > "/dev/stderr"
    bad=1; next
  }
  $2 >= $3 || $4 < $2 || $5 > $3 || $4 >= $5 {
    printf "ERROR: invalid active-HOR/CDR boundaries at line %d.\n", NR > "/dev/stderr"
    bad=1; next
  }
  { print $1, $2, $3, $4, $5 }
  END { exit bad }
' "$activehor" | sort -k1,1 -k2,2n > "$active_clean" ||
  die "Invalid active-HOR/CDR annotation: $activehor"
[[ -s "$active_clean" ]] || die "Active-HOR/CDR annotation contains no intervals."

manifest_clean="${workdir}/fiberseq_manifest.validated.tsv"
awk '
  BEGIN { FS=OFS="\t" }
  /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
  NF != 2 {
    printf "ERROR: manifest line %d has %d fields; expected 2.\n", NR, NF > "/dev/stderr"
    bad=1; next
  }
  { print $1, $2 }
  END { exit bad }
' "$fa_file" > "$manifest_clean" || die "Invalid Fiber-seq manifest: $fa_file"
[[ -s "$manifest_clean" ]] || die "Fiber-seq manifest contains no samples."

duplicate_samples="$(cut -f1 "$manifest_clean" | sort | uniq -d)"
[[ -z "$duplicate_samples" ]] || die $'Duplicate sample IDs found:\n'"$duplicate_samples"

while IFS=$'\t' read -r sample cov; do
  [[ "$sample" =~ ^[A-Za-z0-9._-]+$ ]] ||
    die "Unsafe sample ID '$sample'; use letters, numbers, '.', '_' or '-'."
  [[ -r "$cov" ]] || die "Coverage file for '$sample' not found: $cov"
  [[ -r "${cov}.tbi" || -r "${cov}.csi" ]] ||
    die "Tabix index not found for: $cov"
done < "$manifest_clean"

{
  printf '=========== Call FAS from FIRE coverage ============\n'
  date
  printf 'Fiber-seq manifest: %s\n' "$fa_file"
  printf 'Active-HOR/CDR annotation: %s\n' "$activehor"
  printf 'Output directory: %s\n' "$outdir"
  printf 'Bin size: %s\nSmooth window: %s\nCDR expansion: %s\n' "$winsize" "$smooth_window" "$expansion_size"
  printf 'Outside/inside CDR winsorization: %s/%s\n' "$out_CDR_cutoff" "$in_CDR_cutoff"
  printf 'Threads compatibility value (not used): %s\n' "$threads"
} > "$log_file"

combined_active_cov="${outdir}/FAS_activeHOR_coverage.tsv"
printf 'sample\tchr\tstart\tend\tfire\tlinker\tnucleosome\n' > "$combined_active_cov"
printf '[%s] Extracting active-HOR coverage with tabix...\n' "$(date)" | tee -a "$log_file"

while IFS=$'\t' read -r sample cov; do
  sample_out="${extract_dir}/${sample}.activeHOR.coverage.tsv"
  : > "$sample_out"
  printf '  sample: %s\n' "$sample" | tee -a "$log_file"

  while IFS=$'\t' read -r chr active_start active_end _cdr_start _cdr_end; do
    # Keep the v1.0 region-query behavior. Missing contigs yield no records.
    tabix "$cov" "${chr}:${active_start}-${active_end}" 2>/dev/null || true
  done < "$active_clean" |
    awk -v s="$sample" 'BEGIN{OFS="\t"} NF>=6 {print s,$1,$2,$3,$4,$5,$6}' > "$sample_out"

  cat "$sample_out" >> "$combined_active_cov"
done < "$manifest_clean"

[[ "$(wc -l < "$combined_active_cov")" -gt 1 ]] ||
  die "No active-HOR coverage records were extracted from any sample."

fas_out="${outdir}/FAS_results.tsv"
bin_out="${outdir}/FAS_bin_summary.tsv"
printf '[%s] Running FAS calculation in R...\n' "$(date)" | tee -a "$log_file"

Rscript "$r_code" \
  --coverage "$combined_active_cov" \
  --cen "$active_clean" \
  --out "$fas_out" \
  --bin_out "$bin_out" \
  --plotdir "$plotdir" \
  --bin_size "$winsize" \
  --smooth_k "$smooth_window" \
  --cdr_expansion "$expansion_size" \
  --out_CDR_cutoff "$out_CDR_cutoff" \
  --in_CDR_cutoff "$in_CDR_cutoff" >> "$log_file" 2>&1

printf '[%s] Finished successfully.\n' "$(date)" | tee -a "$log_file"
printf 'FAS results: %s\nBin summary: %s\nQC plots: %s\n' "$fas_out" "$bin_out" "$plotdir"

