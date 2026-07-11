#!/bin/bash
#SBATCH --job-name=call_MARS
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
  bash call_MARS.sh \\
    -ref reference.fa \\
    -cen active_hor_cdr.bed \\
    -o output_dir \\
    -me methylation_files.tsv \\
    -bin 1000 \\
    -n 3 \\
    -dl 5 \\
    -dh 100 \\
    -s 10 \\
    -exp 0 \\
    -out_CDR_cutoff 0.10 \\
    -in_CDR_cutoff 0.05

Required:
  -ref    reference fasta
  -me     four-column TSV: sample_id meth_column depth_column realpath_methylation_file
          methylation file format: chr start end ... methylation_column ... depth_column
          must be bgzip-compressed and tabix-indexed
  -cen    centromere annotation BED:
          chr active_hor_start active_hor_end cdr_start cdr_end

Optional:
  -o      output directory [default: ./output_results]
  -bin    bin size [default: 1000]
  -n      minimum CpG sites per bin [default: 3]
  -dl     minimum CpG depth [default: 5]
  -dh     maximum CpG depth [default: 100]
  -s      smoothing window size [default: 10]
  -exp    symmetric expansion size applied to both CDR boundaries [default: 0]
  -out_CDR_cutoff  winsorization fraction for background outside CDR [default: 0.10]
  -in_CDR_cutoff   winsorization fraction inside CDR [default: 0.05]
EOF
}

ref=""
methfiles=""
cen=""
outdir="output_results"
binsize=1000
cg_cutoff=3
depth_low_cutoff=5
depth_high_cutoff=100
smooth_window=10
expansion_size=0
out_CDR_cutoff=0.10
in_CDR_cutoff=0.05

while [[ $# -gt 0 ]]; do
  case "$1" in
    -ref) ref="$2"; shift 2 ;;
    -me) methfiles="$2"; shift 2 ;;
    -cen) cen="$2"; shift 2 ;;
    -o) outdir="$2"; shift 2 ;;
    -bin) binsize="$2"; shift 2 ;;
    -n) cg_cutoff="$2"; shift 2 ;;
    -dl) depth_low_cutoff="$2"; shift 2 ;;
    -dh) depth_high_cutoff="$2"; shift 2 ;;
    -s) smooth_window="$2"; shift 2 ;;
    -exp) expansion_size="$2"; shift 2 ;;
    -out_CDR_cutoff) out_CDR_cutoff="$2"; shift 2 ;;
    -in_CDR_cutoff) in_CDR_cutoff="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done


## ======================================================================================================================================================================

step_0_check_required_inputs () {
[[ -z "$ref" || -z "$methfiles" || -z "$cen" ]] && { usage; exit 1; }

## check ref file
[[ -f "$ref" ]] || { echo "ERROR: reference not found: $ref"; exit 1; }
[[ "$ref" =~ \.(fa|fasta|fna)(\.gz)?$ ]] || { echo "ERROR: reference should be fa/fasta/fna"; exit 1; }

if [[ ! -f "${ref}.fai" ]]; then
  echo "ERROR: reference index not found: ${ref}.fai"
  echo "Please run: samtools faidx $ref"
  exit 1
fi

## check cen file
[[ -f "$cen" ]] || { echo "ERROR: centromere annotation not found: $cen"; exit 1; }
awk 'NF < 5 {print "ERROR: CEN file must have at least 5 columns at line " NR; exit 1}' "$cen"

## check methfiles
[[ -f "$methfiles" ]] || { echo "ERROR: methylation file list not found: $methfiles"; exit 1; }
awk 'NF < 2 {print "ERROR: methylation file list must have at least two columns at line " NR; exit 1}' "$methfiles"

dup_ids=$(awk '{print $1}' "$methfiles" | sort | uniq -d)
if [[ -n "$dup_ids" ]]; then
	echo "ERROR: duplicated sample IDs: $dup_ids"
	exit 1
fi

while read -r sampleid meth_column depth_column methpath; do
  [[ -f "$methpath" ]] || { echo "ERROR: methylation file not found: $methpath"; exit 1; }
  [[ -f "${methpath}.tbi" || -f "${methpath}.csi" ]] || {
    echo "ERROR: tabix index not found for $methpath"
    exit 1
  }
done < "$methfiles"

## check r code
rcode="$(realpath ./r_code/call_mars.r)"
[[ -f "$rcode" ]] || { echo "ERROR: R code not found: $rcode"; exit 1; }


## check tools
command -v tabix >/dev/null || { echo "ERROR: tabix not found"; exit 1; }
command -v bedtools >/dev/null || { echo "ERROR: bedtools not found"; exit 1; }
command -v Rscript >/dev/null || { echo "ERROR: Rscript not found"; exit 1; }


# Prepare dirs
mkdir -p "$outdir"
outdir=$(realpath "$outdir")
workdir="${outdir}/intermediate_files"
plotdir="${outdir}/MARS_QC_plots"
mkdir -p "$workdir" "$plotdir"

# log file
log="${outdir}/call_MARS.log"

{
  echo "=========== Call MARS ============"
  date
  echo "Reference: $ref"
  echo "Methylation list: $methfiles"
  echo "Centromere annotation: $cen"
  echo "Output directory: $outdir"
  echo "Bin size: $binsize"
  echo "Minimum CpG per bin: $cg_cutoff"
  echo "Depth cutoff: $depth_low_cutoff - $depth_high_cutoff"
  echo "Smooth window: $smooth_window"
  echo "Background trim outside CDR: $out_CDR_cutoff"
  echo "CDR trim inside CDR: $in_CDR_cutoff"
} > "$log"

}

## ======================================================================================================================================================================
step_1_shell_prepare () {

# prepare files
new_methfiles="$workdir/$(basename ${methfiles})"
cp "$methfiles" "$new_methfiles"

modified_cen="$workdir/$(basename ${cen})"
# awk -v expbp="${expansion_size}" '{OFS="\t"; print $1,$2,$3,$4 - expbp,$5 + expbp}' "$cen" |sortBed -i  > $modified_cen

awk -v expand_bp="${expansion_size}" '
{
  OFS = "\t"

  active_start = $2
  active_end   = $3

  cdr_start = $4 - expand_bp
  cdr_end   = $5 + expand_bp

  # positive expansion: do not exceed active HOR boundaries
  if (cdr_start < active_start) cdr_start = active_start
  if (cdr_end > active_end)     cdr_end = active_end

  # also protect against negative coordinates
  if (cdr_start < 0) cdr_start = 0

  # negative expansion: avoid invalid CDR
  if (cdr_end <= cdr_start) {
    printf("ERROR: invalid CDR after boundary adjustment: %s:%d-%d, original CDR: %d-%d, expand_bp: %d\n",
           $1, cdr_start, cdr_end, $4, $5, expand_bp) > "/dev/stderr"
    exit 1
  }

  print $1, active_start, active_end, cdr_start, cdr_end
}
' "$cen" | sortBed -i > "$modified_cen"


# active hor bin windows
binned_active_hor="${workdir}/active_hor_bin${binsize}.bed"
bedtools makewindows -b <(awk '{OFS="\t";print $1,$2,$3}' $modified_cen ) -w "$binsize" | sortBed -i > "${binned_active_hor}"

# extract methylation
combine_binned_meth="${workdir}/combined_active_HOR_bin${binsize}_meth.tsv"

echo -e "chr\tstart\tend\tnCG\tdepth\tmeth\tsample" > "${combine_binned_meth}"

while read -r sampleid meth_column depth_column methpath; do
  [[ -z "${sampleid:-}" ]] && continue
  [[ "$sampleid" =~ ^# ]] && continue
  methpath=$(realpath "$methpath")
  echo -e "\n----- Processing sample: $sampleid" | tee -a "$log"
  echo "  methylation column: $meth_column" | tee -a "$log"
  echo "  depth column: $depth_column" | tee -a "$log"
  echo "  file: $methpath" | tee -a "$log"

  perCG_tmp="${workdir}/${sampleid}_active_hor_perCGDepthFiltered_${depth_low_cutoff}_${depth_high_cutoff}.bed.gz"
  binned_tem="${workdir}/${sampleid}_active_hor_bin${binsize}_meth.bed"

  tabix -R <(awk '{OFS="\t";print $1,$2,$3}' $modified_cen) "$methpath" | awk -v mc="$meth_column" -v dc="$depth_column" -v dl="$depth_low_cutoff" -v dh="$depth_high_cutoff" 'BEGIN{OFS="\t"} $dc >= dl && $dc <= dh {print $1, $2, $3, $mc, $dc }' | sortBed -i |bgzip -f > "$perCG_tmp"

  samscale=$(zcat $perCG_tmp |cut -f 4 | shuf -n 1000 | awk ' $1 != "NA" {sum += $1; n++ } END {if (n == 0) print "NA";else print sum/n;}')

  echo -e "  sampled methylation mean by shuf 1000 rows: ${samscale}\n" | tee -a "$log"

  if [[ "$samscale" == "NA" ]]; then
	  echo "ERROR: no valid methylation values found for $sampleid" | tee -a "$log"
	  exit 1
  fi

  bedtools map -sorted -a "${binned_active_hor}" -b "$perCG_tmp" -c 4,5,4 -o mean,mean,count -null NA | awk -v scale="${samscale}" -v cg="$cg_cutoff" -v sid="$sampleid" 'BEGIN{OFS="\t"} $6 != "NA" && $6 >= cg {if (scale <= 1) $4=$4*100 ;printf("%s\t%s\t%s\t%s\t%.1f\t%.4f\t%s\n", $1,$2,$3,$6,$5,$4,sid)}' > "$binned_tem"
  cat "$binned_tem" >> "${combine_binned_meth}"
done < "$new_methfiles"

bgzip -f "${combine_binned_meth}"

echo -e "\n\n---------- step 1 finished, combined binned methylation written to: ${combine_binned_meth}.gz" | tee -a "$log"
}

## ======================================================================================================================================================================

step_2_run_r () {
mars_out="${outdir}/MARS_results.tsv"

Rscript "$rcode" \
  --meth "${combine_binned_meth}.gz" \
  --cen "${modified_cen}" \
  --out "$mars_out" \
  --plotdir "$plotdir" \
  --smooth_k "$smooth_window" \
  --out_CDR_cutoff "$out_CDR_cutoff" \
  --in_CDR_cutoff "$in_CDR_cutoff" \
  --bin_size "$binsize" \
  >> "$log" 2>&1
}

# main
step_0_check_required_inputs
step_1_shell_prepare
step_2_run_r

{
  echo "All done."
  date
} | tee -a "$log"



