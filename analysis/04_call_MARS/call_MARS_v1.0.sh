#!/usr/bin/env bash
# Calculate MARS from methylation data using a fixed active-HOR/CDR annotation.
# The embedded Slurm settings are ignored when the script is run with bash.
#SBATCH --job-name=call_MARS
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

set -Eeuo pipefail

readonly PROGRAM_NAME="$(basename "$0")"
readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"

usage() {
  cat <<EOF
Usage:
  ${PROGRAM_NAME} --centromere FILE --methylation-list FILE [options]

Required arguments:
  --centromere, -cen FILE         Five-column active-HOR/CDR annotation
  --methylation-list, -me FILE    Four-column sample manifest

Analysis options (defaults reproduce v1.0):
  --output-dir, -o DIR            Output directory [output_results]
  --bin-size, -bin INT            Bin size in bp [1000]
  --min-cpg, -n INT               Minimum CpG sites per bin [3]
  --min-depth, -dl NUMBER         Minimum per-CpG depth [5]
  --max-depth, -dh NUMBER         Maximum per-CpG depth [100]
  --smooth-window, -s INT         Smoothing window in bins [10]
  --cdr-expansion, -exp INT       Symmetric CDR expansion in bp [0]
  --outside-cdr-cutoff, -out_CDR_cutoff FRACTION
                                      Outside-CDR winsorization fraction [0.10]
  --inside-cdr-cutoff, -in_CDR_cutoff FRACTION
                                      Inside-CDR winsorization fraction [0.05]
  -h, --help                      Show this message and exit

Centromere annotation (tab-delimited):
  chr  active_hor_start  active_hor_end  CDR_start  CDR_end

Sample manifest (tab-delimited):
  sample_id  methylation_column  depth_column  bgzip_methylation_file

Comment lines beginning with '#' and blank lines are ignored. Column numbers
are one-based. Relative methylation paths are resolved from the launch directory.

Example:
  bash ${PROGRAM_NAME} \\
    --centromere example_data/activeHOR_CDR.tsv \\
    --methylation-list example_data/meth_files.tsv \\
    --output-dir example_output
EOF
}

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
log_msg() { printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*" | tee -a "$log_file" >&2; }
require_value() { [[ $# -ge 2 && -n "${2:-}" ]] || die "Option '$1' requires a value."; }
require_command() { command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"; }
is_positive_integer() { [[ "$1" =~ ^[1-9][0-9]*$ ]]; }
is_integer() { [[ "$1" =~ ^-?[0-9]+$ ]]; }
is_nonnegative_number() { [[ "$1" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]]; }
is_fraction() { is_nonnegative_number "$1" && awk -v x="$1" 'BEGIN { exit !(x < 0.5) }'; }

methylation_list=""
centromere=""
output_dir="output_results"
bin_size=1000
min_cpg=3
min_depth=5
max_depth=100
smooth_window=10
cdr_expansion=0
outside_cdr_cutoff=0.10
inside_cdr_cutoff=0.05

while (($#)); do
  case "$1" in
    --methylation-list|-me) require_value "$@"; methylation_list="$2"; shift 2 ;;
    --centromere|-cen) require_value "$@"; centromere="$2"; shift 2 ;;
    --output-dir|-o) require_value "$@"; output_dir="$2"; shift 2 ;;
    --bin-size|-bin) require_value "$@"; bin_size="$2"; shift 2 ;;
    --min-cpg|-n) require_value "$@"; min_cpg="$2"; shift 2 ;;
    --min-depth|-dl) require_value "$@"; min_depth="$2"; shift 2 ;;
    --max-depth|-dh) require_value "$@"; max_depth="$2"; shift 2 ;;
    --smooth-window|-s) require_value "$@"; smooth_window="$2"; shift 2 ;;
    --cdr-expansion|-exp) require_value "$@"; cdr_expansion="$2"; shift 2 ;;
    --outside-cdr-cutoff|-out_CDR_cutoff) require_value "$@"; outside_cdr_cutoff="$2"; shift 2 ;;
    --inside-cdr-cutoff|-in_CDR_cutoff) require_value "$@"; inside_cdr_cutoff="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    --) shift; break ;;
    *) die "Unknown option: $1 (run '${PROGRAM_NAME} --help')" ;;
  esac
done
(($# == 0)) || die "Unexpected positional argument(s): $*"

[[ -n "$methylation_list" ]] || die "--methylation-list is required."
[[ -n "$centromere" ]] || die "--centromere is required."
[[ -r "$methylation_list" ]] || die "Methylation manifest not found or unreadable: $methylation_list"
[[ -r "$centromere" ]] || die "Centromere annotation not found or unreadable: $centromere"

is_positive_integer "$bin_size" || die "--bin-size must be a positive integer."
is_positive_integer "$min_cpg" || die "--min-cpg must be a positive integer."
is_nonnegative_number "$min_depth" || die "--min-depth must be non-negative."
is_nonnegative_number "$max_depth" || die "--max-depth must be non-negative."
awk -v lo="$min_depth" -v hi="$max_depth" 'BEGIN { exit !(lo <= hi) }' || die "--min-depth cannot exceed --max-depth."
is_positive_integer "$smooth_window" || die "--smooth-window must be a positive integer."
is_integer "$cdr_expansion" || die "--cdr-expansion must be an integer (negative values contract the CDR)."
is_fraction "$outside_cdr_cutoff" || die "--outside-cdr-cutoff must be in [0, 0.5)."
is_fraction "$inside_cdr_cutoff" || die "--inside-cdr-cutoff must be in [0, 0.5)."

for tool in awk bedtools bgzip cut Rscript shuf sort tabix zcat; do require_command "$tool"; done

r_code="${SCRIPT_DIR}/r_code/call_mars.r"
[[ -r "$r_code" ]] || die "R program not found: $r_code"

mkdir -p -- "$output_dir"
output_dir="$(cd -- "$output_dir" && pwd -P)"
work_dir="${output_dir}/intermediate_files"
plot_dir="${output_dir}/MARS_QC_plots"
mkdir -p -- "$work_dir" "$plot_dir"
log_file="${output_dir}/call_MARS.log"

clean_manifest="${work_dir}/methylation_manifest.validated.tsv"
awk '
  BEGIN { FS=OFS="\t" }
  /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
  NF != 4 {
    printf "ERROR: manifest line %d has %d fields; expected 4.\n", NR, NF > "/dev/stderr"
    bad=1; next
  }
  { print }
  END { exit bad }
' "$methylation_list" > "$clean_manifest" || die "Invalid methylation manifest: $methylation_list"
[[ -s "$clean_manifest" ]] || die "Methylation manifest contains no samples."

duplicates="$(cut -f1 "$clean_manifest" | sort | uniq -d)"
[[ -z "$duplicates" ]] || die $'Duplicate sample IDs found:\n'"$duplicates"
while IFS=$'\t' read -r sample_id meth_column depth_column meth_path; do
  [[ "$sample_id" =~ ^[A-Za-z0-9._-]+$ ]] || die "Unsafe sample ID '$sample_id'; allowed: letters, numbers, '.', '_' and '-'."
  is_positive_integer "$meth_column" || die "Invalid methylation column for '$sample_id': $meth_column"
  is_positive_integer "$depth_column" || die "Invalid depth column for '$sample_id': $depth_column"
  [[ "$meth_column" != "$depth_column" ]] || die "Methylation and depth columns are identical for '$sample_id'."
  [[ -r "$meth_path" ]] || die "Methylation file for '$sample_id' not found: $meth_path"
  [[ -r "${meth_path}.tbi" || -r "${meth_path}.csi" ]] || die "Tabix index missing for: $meth_path"
done < "$clean_manifest"

validated_cen="${work_dir}/active_HOR_CDR.validated.tsv"
awk '
  BEGIN { FS=OFS="\t" }
  /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
  NF < 5 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $4 !~ /^[0-9]+$/ || $5 !~ /^[0-9]+$/ {
    printf "ERROR: invalid five-column annotation at line %d.\n", NR > "/dev/stderr"; bad=1; next
  }
  $2 >= $3 || $4 < $2 || $5 > $3 || $4 >= $5 {
    printf "ERROR: invalid active-HOR/CDR boundaries at line %d.\n", NR > "/dev/stderr"; bad=1; next
  }
  { print $1, $2, $3, $4, $5 }
  END { exit bad }
' "$centromere" | sort -k1,1 -k2,2n > "$validated_cen" || die "Invalid centromere annotation: $centromere"
[[ -s "$validated_cen" ]] || die "Centromere annotation contains no intervals."

modified_cen="${work_dir}/active_HOR_CDR.expanded_${cdr_expansion}bp.tsv"
awk -v expand_bp="$cdr_expansion" '
  BEGIN { OFS="\t" }
  {
    cdr_start=$4-expand_bp; cdr_end=$5+expand_bp
    if (cdr_start<$2) cdr_start=$2
    if (cdr_end>$3) cdr_end=$3
    if (cdr_start<0) cdr_start=0
    if (cdr_end<=cdr_start) {
      printf "ERROR: invalid CDR after expansion: %s:%d-%d\n", $1, cdr_start, cdr_end > "/dev/stderr"; exit 1
    }
    print $1, $2, $3, cdr_start, cdr_end
  }
' "$validated_cen" > "$modified_cen"

region_bed="${work_dir}/active_HOR.regions.bed"
cut -f1-3 "$modified_cen" > "$region_bed"
binned_active_hor="${work_dir}/active_HOR.${bin_size}bp_bins.bed"
bedtools makewindows -b "$region_bed" -w "$bin_size" | sort -k1,1 -k2,2n > "$binned_active_hor"

combined_tsv="${work_dir}/combined_active_HOR.${bin_size}bp_methylation.tsv"
printf 'chr\tstart\tend\tnCG\tdepth\tmeth\tsample\n' > "$combined_tsv"

{
  printf '=========== Call MARS ============\n'
  date
  printf 'Manifest: %s\nAnnotation: %s\nOutput: %s\n' "$methylation_list" "$centromere" "$output_dir"
  printf 'Bin size: %s\nMinimum CpG/bin: %s\nDepth: %s-%s\nSmooth window: %s\nCDR expansion: %s\n' "$bin_size" "$min_cpg" "$min_depth" "$max_depth" "$smooth_window" "$cdr_expansion"
  printf 'Outside/inside CDR winsorization: %s/%s\n' "$outside_cdr_cutoff" "$inside_cdr_cutoff"
} > "$log_file"

sample_count="$(wc -l < "$clean_manifest" | tr -d ' ')"
log_msg "Validated ${sample_count} sample(s)."
while IFS=$'\t' read -r sample_id meth_column depth_column meth_path; do
  log_msg "Processing ${sample_id}"
  per_cpg="${work_dir}/${sample_id}.active_HOR.depth_${min_depth}-${max_depth}.bed.gz"
  binned_sample="${work_dir}/${sample_id}.active_HOR.${bin_size}bp_methylation.tsv"

  tabix -R "$region_bed" "$meth_path" |
    awk -v mc="$meth_column" -v dc="$depth_column" -v dl="$min_depth" -v dh="$max_depth" '
      BEGIN { OFS="\t" }
      NF < mc || NF < dc { print "ERROR: methylation record has too few columns" > "/dev/stderr"; exit 2 }
      $dc >= dl && $dc <= dh { print $1, $2, $3, $mc, $dc }
    ' |
    sort -k1,1 -k2,2n | bgzip -c > "$per_cpg"

  sample_scale="$(zcat "$per_cpg" | cut -f4 | shuf -n 1000 | awk '$1 != "NA" && $1 != "." && $1 != "" {sum+=$1; n++} END {if (!n) print "NA"; else print sum/n}')"
  [[ "$sample_scale" != "NA" ]] || die "No valid methylation values found for '$sample_id' after depth filtering."
  log_msg "${sample_id}: mean of up to 1,000 randomly sampled methylation values = ${sample_scale}"

  bedtools map -sorted -a "$binned_active_hor" -b "$per_cpg" -c 4,5,4 -o mean,mean,count -null NA |
    awk -v scale="$sample_scale" -v minimum="$min_cpg" -v sid="$sample_id" '
      BEGIN { OFS="\t" }
      $6 != "NA" && $6 >= minimum {
        if (scale <= 1) $4=$4*100
        printf "%s\t%s\t%s\t%s\t%.1f\t%.4f\t%s\n", $1,$2,$3,$6,$5,$4,sid
      }
    ' > "$binned_sample"
  cat "$binned_sample" >> "$combined_tsv"
done < "$clean_manifest"

bgzip -f "$combined_tsv"
mars_output="${output_dir}/MARS_results.tsv"
log_msg "Running MARS calculation."
Rscript "$r_code" \
  --meth "${combined_tsv}.gz" --cen "$modified_cen" --out "$mars_output" --plotdir "$plot_dir" \
  --smooth_k "$smooth_window" --out_CDR_cutoff "$outside_cdr_cutoff" \
  --in_CDR_cutoff "$inside_cdr_cutoff" --bin_size "$bin_size" >> "$log_file" 2>&1

log_msg "Finished successfully: ${mars_output}"

