#!/usr/bin/env bash
# Call candidate centromeric dip regions (CDRs) from binned methylation data.
#
# This script can be submitted with sbatch or run directly with bash.
#SBATCH --job-name=call_CDR
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
  ${PROGRAM_NAME} --active-hor FILE --methylation-list FILE [options]

Required arguments:
  --active-hor, -cen FILE       BED file containing active HOR regions
  --methylation-list, -me FILE  Four-column, tab-delimited sample manifest

Output:
  --output-dir, -o DIR          Output directory (default: CDR_output)

CDR parameters:
  --bin-size INT                Window size in bp (default: 5000)
  --bin-gap INT                 Maximum gap between bins in bp (default: 50000)
  --min-bins INT                Minimum number of bins per candidate (default: 5)
  --cutoff-frac FLOAT           Methylation cutoff fraction (default: 0.15)
  --smooth-k INT                Smoothing parameter (default: 10)
  --depth-cutoff FLOAT          Minimum mean depth (default: 0)
  --min-sites-per-bin INT       Minimum mapped records per bin (default: 3)

Other options:
  --keep-intermediates          Retain per-sample binned methylation tables
  -h, --help                    Show this help message and exit

Manifest format (tab-delimited; comments beginning with '#' are ignored):
  sample_name  methylation_column  depth_column  bgzip_methylation_file

Column numbers are one-based. Methylation files must be coordinate-sorted,
bgzip-compressed, and indexed with tabix. Relative data paths are resolved from
the directory in which this command is launched.

Example:
  bash ${PROGRAM_NAME} \\
    --active-hor example_data/PAN027_Mat_chr9_activeHOR_annotation.tsv \\
    --methylation-list example_data/meth_file.tsv \\
    --output-dir example_output
EOF
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

log() {
    printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*" >&2
}

require_value() {
    [[ $# -ge 2 && -n "${2:-}" ]] || die "Option '$1' requires a value."
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found in PATH: $1"
}

is_positive_integer() {
    [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

is_nonnegative_number() {
    [[ "$1" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]]
}

active_hor=""
output_dir="CDR_output"
methylation_list=""
bin_size=5000
bin_gap=50000
min_bins=5
cutoff_frac=0.15
smooth_k=10
depth_cutoff=0
min_sites_per_bin=3
keep_intermediates=false

while (($# > 0)); do
    case "$1" in
        --active-hor|-cen)
            require_value "$@"; active_hor="$2"; shift 2 ;;
        --output-dir|-o)
            require_value "$@"; output_dir="$2"; shift 2 ;;
        --methylation-list|-me)
            require_value "$@"; methylation_list="$2"; shift 2 ;;
        --bin-size|-bin_size)
            require_value "$@"; bin_size="$2"; shift 2 ;;
        --bin-gap|-bin_gap)
            require_value "$@"; bin_gap="$2"; shift 2 ;;
        --min-bins|-min_bins)
            require_value "$@"; min_bins="$2"; shift 2 ;;
        --cutoff-frac|-cutoff_frac)
            require_value "$@"; cutoff_frac="$2"; shift 2 ;;
        --smooth-k|-smooth_k)
            require_value "$@"; smooth_k="$2"; shift 2 ;;
        --depth-cutoff|-depth_cutoff)
            require_value "$@"; depth_cutoff="$2"; shift 2 ;;
        --min-sites-per-bin)
            require_value "$@"; min_sites_per_bin="$2"; shift 2 ;;
        --keep-intermediates)
            keep_intermediates=true; shift ;;
        -h|--help)
            usage; exit 0 ;;
        --)
            shift; break ;;
        *)
            die "Unknown option: $1 (run '${PROGRAM_NAME} --help' for usage)" ;;
    esac
done

(($# == 0)) || die "Unexpected positional argument(s): $*"
[[ -n "$active_hor" ]] || die "--active-hor is required."
[[ -n "$methylation_list" ]] || die "--methylation-list is required."
[[ -f "$active_hor" ]] || die "Active HOR BED file not found: $active_hor"
[[ -r "$active_hor" ]] || die "Active HOR BED file is not readable: $active_hor"
[[ -f "$methylation_list" ]] || die "Methylation manifest not found: $methylation_list"
[[ -r "$methylation_list" ]] || die "Methylation manifest is not readable: $methylation_list"

is_positive_integer "$bin_size" || die "--bin-size must be a positive integer."
is_positive_integer "$bin_gap" || die "--bin-gap must be a positive integer."
is_positive_integer "$min_bins" || die "--min-bins must be a positive integer."
is_positive_integer "$smooth_k" || die "--smooth-k must be a positive integer."
is_positive_integer "$min_sites_per_bin" || die "--min-sites-per-bin must be a positive integer."
is_nonnegative_number "$cutoff_frac" || die "--cutoff-frac must be a non-negative number."
is_nonnegative_number "$depth_cutoff" || die "--depth-cutoff must be a non-negative number."
awk -v x="$cutoff_frac" 'BEGIN { exit !(x <= 1) }' || die "--cutoff-frac must be between 0 and 1."

for command_name in awk bedtools Rscript sort tabix; do
    require_command "$command_name"
done

r_script="${SCRIPT_DIR}/r_code/call_CDR_from_active_HOR.R"
[[ -f "$r_script" ]] || die "R script not found: $r_script"

mkdir -p -- "$output_dir"
output_dir="$(cd -- "$output_dir" && pwd -P)"
work_dir="$(mktemp -d "${TMPDIR:-/tmp}/call_CDR.XXXXXX")"
cleanup() {
    rm -rf -- "$work_dir"
}
trap cleanup EXIT
trap 'die "Command failed at line ${LINENO}."' ERR

manifest_clean="${work_dir}/manifest.tsv"
awk '
    BEGIN { FS=OFS="\t" }
    /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
    NF != 4 {
        printf "ERROR: manifest line %d has %d fields; expected 4.\n", NR, NF > "/dev/stderr"
        bad=1; next
    }
    { print }
    END { exit bad }
' "$methylation_list" > "$manifest_clean" || die "Invalid methylation manifest: $methylation_list"
[[ -s "$manifest_clean" ]] || die "The methylation manifest contains no samples."

duplicate_samples="$(cut -f1 "$manifest_clean" | sort | uniq -d)"
[[ -z "$duplicate_samples" ]] || die $'Duplicate sample names found:\n'"$duplicate_samples"

while IFS=$'\t' read -r sample meth_col depth_col meth_file; do
    [[ "$sample" =~ ^[A-Za-z0-9._-]+$ ]] || \
        die "Invalid sample name '$sample'; use only letters, numbers, '.', '_' and '-'."
    is_positive_integer "$meth_col" || die "Methylation column for '$sample' is not a positive integer: $meth_col"
    is_positive_integer "$depth_col" || die "Depth column for '$sample' is not a positive integer: $depth_col"
    [[ "$meth_col" != "$depth_col" ]] || die "Methylation and depth columns are identical for '$sample'."
    [[ -f "$meth_file" ]] || die "Methylation file for '$sample' not found: $meth_file"
    [[ -r "$meth_file" ]] || die "Methylation file for '$sample' is not readable: $meth_file"
    [[ -f "${meth_file}.tbi" || -f "${meth_file}.csi" ]] || \
        die "Tabix index not found for '$meth_file' (expected .tbi or .csi)."
done < "$manifest_clean"

region_bed="${work_dir}/active_HOR.regions.bed"
awk '
    BEGIN { FS=OFS="\t" }
    /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
    NF < 3 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $2 >= $3 {
        printf "ERROR: invalid BED interval at active-HOR line %d.\n", NR > "/dev/stderr"
        bad=1; next
    }
    { print $1, $2, $3 }
    END { exit bad }
' "$active_hor" | sort -k1,1 -k2,2n > "$region_bed" || die "Invalid active HOR BED file: $active_hor"
[[ -s "$region_bed" ]] || die "The active HOR BED file contains no intervals."

binned_bed="${work_dir}/active_HOR.regions_bin${bin_size}.bed"
bedtools makewindows -b "$region_bed" -w "$bin_size" | \
    sort -k1,1 -k2,2n > "$binned_bed"

sample_count="$(wc -l < "$manifest_clean" | tr -d ' ')"
log "Validated ${sample_count} sample(s); results will be written to ${output_dir}"

while IFS=$'\t' read -r sample meth_col depth_col meth_file; do
    log "Processing ${sample}"
    extracted="${work_dir}/${sample}.active_HOR.meth.tsv"
    mapped="${work_dir}/${sample}.mapped.tsv"

    printf 'chr\tstart\tend\tmeth\tdepth\n' > "$extracted"
    tabix -R "$region_bed" "$meth_file" | \
        awk -v mc="$meth_col" -v dc="$depth_col" '
            BEGIN { OFS="\t" }
            NF < mc || NF < dc {
                printf "ERROR: input record has fewer than the requested columns.\n" > "/dev/stderr"
                exit 2
            }
            { print $1, $2, $3, $mc, $dc }
        ' | \
        sort -k1,1 -k2,2n | \
        bedtools map -a "$binned_bed" -b - -c 4,5,5 -o mean,mean,count > "$mapped"

    awk -v minimum="$min_sites_per_bin" '
        BEGIN { OFS="\t" }
        $6 >= minimum { print $1, $2, $3, $4, $5 }
    ' "$mapped" >> "$extracted"

    Rscript "$r_script" \
        --sample "$sample" \
        --meth_file "$extracted" \
        --active_hor "$active_hor" \
        --out_dir "$output_dir" \
        --bin_gap "$bin_gap" \
        --min_bins "$min_bins" \
        --cutoff_frac "$cutoff_frac" \
        --smooth_k "$smooth_k" \
        --depth_cutoff "$depth_cutoff"

    if [[ "$keep_intermediates" == true ]]; then
        cp -- "$extracted" "${output_dir}/${sample}.active_HOR.meth.tsv"
    fi
done < "$manifest_clean"

log "Finished successfully."



