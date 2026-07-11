#!/bin/bash
#SBATCH --job-name=call_CDR
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

set -euo pipefail

usage() {
cat <<EOF
Usage:
  bash call_CDR_from_active_HOR.sh \\
    -cen active_hor.bed \\
    -o output_dir \\
    -me methylation_files.tsv \\
    -bin_size 5000 \\
    -bin_gap 50000 \\
    -min_bins 3 \\
    -cutoff_frac 0.15 \\
    -smooth_k 10

Input methylation_files.tsv: 4 columns, no header
  sample_name    meth_col    depth_col    bgzip_tabix_meth_file

Example:
  PAN010_Mat    meth    depth    /path/PAN010_Mat.meth.bed.gz

Required:
  methylation files must be bgzip compressed and tabix indexed.
  sample_name must be unique.
EOF
}

# -------------------------
# default parameters
# -------------------------
cen=""
outdir="CDR_output"
meth_list=""
bin_size=5000
bin_gap=50000
min_bins=5
cutoff_frac=0.15
smooth_k=10
depth_cutoff=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -cen) cen="$2"; shift 2 ;;
        -o) outdir="$2"; shift 2 ;;
        -me) meth_list="$2"; shift 2 ;;
	-bin_size) bin_size="$2"; shift 2 ;;
        -bin_gap) bin_gap="$2"; shift 2 ;;
        -min_bins) min_bins="$2"; shift 2 ;;
        -cutoff_frac) cutoff_frac="$2"; shift 2 ;;
        -smooth_k) smooth_k="$2"; shift 2 ;;
        -depth_cutoff) depth_cutoff="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1"; usage; exit 1 ;;
    esac
done

[[ -n "$cen" && -n "$meth_list" ]] || { usage; exit 1; }
[[ -f "$cen" ]] || { echo "ERROR: active HOR BED not found: $cen"; exit 1; }
[[ -f "$meth_list" ]] || { echo "ERROR: methylation list not found: $meth_list"; exit 1; }

mkdir -p "$outdir"

# -------------------------
# check duplicated sample names
# -------------------------
dup=$(awk '{print $1}' "$meth_list" | sort | uniq -d)
if [[ -n "$dup" ]]; then
    echo "ERROR: duplicated sample names found:"
    echo "$dup"
    exit 1
fi

# -------------------------
# check bgzip/tabix index
# -------------------------
while read -r sample meth_col depth_col meth_file; do
    [[ -f "$meth_file" ]] || { echo "ERROR: file not found: $meth_file"; exit 1; }
    [[ -f "${meth_file}.tbi" || -f "${meth_file}.csi" ]] || {
        echo "ERROR: tabix index not found for $meth_file"
        echo "Please run: tabix -p bed $meth_file"
        exit 1
    }
done < "$meth_list"

# -------------------------
# extract active HOR regions and call CDR sample by sample
# -------------------------
r_script=$(realpath call_CDR_from_active_HOR.R)

while read -r sample meth_col depth_col meth_file; do
    echo "## Processing $sample"

    region_bed="active_HOR.regions.bed"
    awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' "$cen" > "$region_bed"
    region_bed_binned="active_HOR.regions_bin${bin_size}.bed"
    
    bedtools makewindows -b $region_bed -w "$bin_size" | sortBed -i > "${region_bed_binned}"

    extracted="${sample}.active_HOR.meth.tsv"
    echo -e "chr\tstart\tend\tmeth\tdepth" > $extracted
    tabix -R "$region_bed" "$meth_file" | awk -v mc="$meth_col" -v dc="$depth_col" ' {OFS="\t"; print $1, $2, $3, $mc, $dc }' |bedtools map -a $region_bed_binned -b - -c 4,5,5 -o mean,mean,count |awk '$NF >=3 {OFS="\t";print $1,$2,$3,$4,$5}'  >> "$extracted"

    Rscript "$r_script" \
        --sample "$sample" \
        --meth_file "$extracted" \
        --active_hor "$cen" \
        --out_dir "$outdir" \
        --bin_gap "$bin_gap" \
        --min_bins "$min_bins" \
        --cutoff_frac "$cutoff_frac" \
        --smooth_k "$smooth_k" \
        --depth_cutoff "$depth_cutoff"

done < "$meth_list"

echo "Done."


