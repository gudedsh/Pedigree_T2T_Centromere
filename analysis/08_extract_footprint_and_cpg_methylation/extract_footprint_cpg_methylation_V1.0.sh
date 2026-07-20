#!/usr/bin/env bash
# Extract inferred nucleosome footprints and CpG methylation from a tagged BAM.
#SBATCH --job-name=footprint_cpg
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=48:00:00

set -Eeuo pipefail

readonly PROGRAM_NAME="$(basename "$0")"
readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"

usage() {
  cat <<'EOF'
Usage:
  extract_footprint_cpg_methylation.sh \
    --bam FILE \
    --regions FILE \
    --reference FILE \
    --sample NAME [options]

Required:
  --bam FILE                 Aligned BAM containing Fiber-seq tags
  --regions FILE             Five-column active-HOR/CDR annotation
  --reference FILE           Reference FASTA used for BAM alignment
  --sample NAME              Sample ID, e.g. PAN027Pat_NPC

Options:
  --assembly NAME            Assembly label [sample name before final underscore]
  --cell NAME                Cell label [sample name after final underscore]
  --output-dir DIR           Output directory [footprint_cpg_output]
  --threads INT              Worker threads [SLURM_CPUS_PER_TASK or 8]
  --mapq INT                 Minimum read mapping quality [1]
  --keep-intermediates       Retain all intermediate BED and BAM files
  -h, --help                 Show this help message

Example:
  bash extract_footprint_cpg_methylation.sh \
    --bam example_data/PAN027Pat_chr9_activeHOR_NPC_aligned.bam \
    --regions example_data/activeHOR_CDR.tsv \
    --reference example_data/PAN027Pat_chr9.fa.gz \
    --sample PAN027Pat_NPC \
    --output-dir example_output
EOF
}

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
log() { printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*" >&2; }
require_value() { [[ $# -ge 2 && -n "${2:-}" ]] || die "Option '$1' requires a value."; }
require_file() { [[ -s "$1" ]] || die "File not found or empty: $1"; }
require_command() { command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"; }
is_nonnegative_integer() { [[ "$1" =~ ^[0-9]+$ ]]; }
is_positive_integer() { [[ "$1" =~ ^[1-9][0-9]*$ ]]; }

bam=""
regions=""
reference=""
sample=""
assembly=""
cell=""
output_dir="footprint_cpg_output"
threads="${SLURM_CPUS_PER_TASK:-8}"
mapq=1
keep_intermediates=false

while (($#)); do
  case "$1" in
    --bam) require_value "$@"; bam="$2"; shift 2 ;;
    --regions|--centromere) require_value "$@"; regions="$2"; shift 2 ;;
    --reference) require_value "$@"; reference="$2"; shift 2 ;;
    --sample) require_value "$@"; sample="$2"; shift 2 ;;
    --assembly) require_value "$@"; assembly="$2"; shift 2 ;;
    --cell) require_value "$@"; cell="$2"; shift 2 ;;
    --output-dir|-o) require_value "$@"; output_dir="$2"; shift 2 ;;
    --threads|-t) require_value "$@"; threads="$2"; shift 2 ;;
    --mapq) require_value "$@"; mapq="$2"; shift 2 ;;
    --keep-intermediates) keep_intermediates=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1 (run '${PROGRAM_NAME} --help')" ;;
  esac
done

[[ -n "$bam" ]] || die "--bam is required."
[[ -n "$regions" ]] || die "--regions is required."
[[ -n "$reference" ]] || die "--reference is required."
[[ -n "$sample" ]] || die "--sample is required."
require_file "$bam"
require_file "$regions"
require_file "$reference"
[[ -s "${bam}.bai" || -s "${bam%.bam}.bai" || -s "${bam}.csi" ]] ||
  die "BAM index not found for: $bam"
[[ "$sample" =~ ^[A-Za-z0-9._-]+$ ]] || die "Unsafe sample ID: $sample"
is_positive_integer "$threads" || die "--threads must be a positive integer."
is_nonnegative_integer "$mapq" || die "--mapq must be a non-negative integer."

for tool in awk bedtools bgzip ft gzip samtools sort; do require_command "$tool"; done

if [[ -z "$assembly" ]]; then
  [[ "$sample" == *_* ]] || die "Cannot infer assembly from '$sample'; provide --assembly."
  assembly="${sample%_*}"
fi
if [[ -z "$cell" ]]; then
  [[ "$sample" == *_* ]] || die "Cannot infer cell from '$sample'; provide --cell."
  cell="${sample##*_}"
fi

mkdir -p -- "$output_dir"
output_dir="$(cd -- "$output_dir" && pwd -P)"
work_dir="${output_dir}/intermediate_files"
reference_cache_dir="${work_dir}/reference_cache"
sort_tmp="${work_dir}/sort_tmp"
mkdir -p -- "$work_dir" "$reference_cache_dir" "$sort_tmp"
log_file="${output_dir}/extract_footprint_cpg_methylation.log"

sort_bed() {
  LC_ALL=C sort -k1,1 -k2,2n --parallel="$threads" -S 4G -T "$sort_tmp"
}

bed12_blocks_to_bed6() {
  local input_bed=$1
  local output_bed=$2
  require_file "$input_bed"

  awk '
    BEGIN { OFS="\t" }
    NF < 12 { printf "ERROR: malformed BED12 record at line %d.\n", NR > "/dev/stderr"; bad=1; next }
    {
      split($11, block_sizes, ",")
      split($12, block_starts, ",")
      for (i=1; i<=$10; i++) {
        start=$2+block_starts[i]
        end=start+block_sizes[i]
        if (end<=start) end=start+1
        print $1,start,end,$4,".",$6
      }
    }
    END { exit bad }
  ' "$input_bed" | sort_bed > "$output_bed"
}

# Cache an uncompressed FASTA because some BEDTools builds segfault on .fa.gz.
if [[ "$reference" == *.gz ]]; then
  reference_for_tools="${reference_cache_dir}/$(basename "${reference%.gz}")"
  if [[ ! -s "$reference_for_tools" || "$reference" -nt "$reference_for_tools" ]]; then
    log "Decompressing reference for BEDTools and fibertools."
    gzip -dc -- "$reference" > "${reference_for_tools}.tmp"
    require_file "${reference_for_tools}.tmp"
    mv -- "${reference_for_tools}.tmp" "$reference_for_tools"
  fi
else
  reference_for_tools="$reference"
fi
if [[ ! -s "${reference_for_tools}.fai" || "$reference_for_tools" -nt "${reference_for_tools}.fai" ]]; then
  samtools faidx "$reference_for_tools"
fi
require_file "${reference_for_tools}.fai"

active_bed="${work_dir}/${sample}.active_HOR.bed"
cdr_bed="${work_dir}/${sample}.CDR.bed"
active_bam="${work_dir}/${sample}.active_HOR.bam"
cpg_bed12="${work_dir}/${sample}.ft_extract_cpg.bed"
nuc_bed12="${work_dir}/${sample}.ft_extract_nuc.bed"
single_cpg_bed="${work_dir}/${sample}.single_cpg.bed"
single_nuc_bed="${work_dir}/${sample}.single_nucleosome.bed"
nuc_cpg_content="${work_dir}/${sample}.nucleosome_reference_CpG.bed"
methylated_counts="${work_dir}/${sample}.nucleosome_methylated_CpG_counts.bed"
final_tsv="${output_dir}/${sample}.nucleosome_CpG_methylation.tsv"
final_gz="${final_tsv}.gz"

{
  printf '=========== Extract nucleosome footprints and CpG methylation ============\n'
  date
  printf 'Sample: %s\nAssembly: %s\nCell: %s\n' "$sample" "$assembly" "$cell"
  printf 'BAM: %s\nRegions: %s\nReference: %s\n' "$bam" "$regions" "$reference"
  printf 'Reference used by tools: %s\nOutput: %s\nThreads: %s\nMAPQ: %s\n'     "$reference_for_tools" "$output_dir" "$threads" "$mapq"
} > "$log_file"

log "Preparing active-HOR and CDR intervals."
: > "$active_bed"
: > "$cdr_bed"
awk '
  BEGIN { OFS="\t" }
  /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
  NF < 5 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $4 !~ /^[0-9]+$/ || $5 !~ /^[0-9]+$/ {
    printf "ERROR: invalid five-column annotation at line %d.\n", NR > "/dev/stderr"; bad=1; next
  }
  $2 >= $3 || $4 < $2 || $5 > $3 || $4 >= $5 {
    printf "ERROR: invalid active-HOR/CDR boundaries at line %d.\n", NR > "/dev/stderr"; bad=1; next
  }
  { print $1,$2,$3 > active; print $1,$4,$5 > cdr }
  END { exit bad }
' active="$active_bed" cdr="$cdr_bed" "$regions"
require_file "$active_bed"
require_file "$cdr_bed"
sort_bed < "$active_bed" > "${active_bed}.sorted"
sort_bed < "$cdr_bed" > "${cdr_bed}.sorted"
mv -- "${active_bed}.sorted" "$active_bed"
mv -- "${cdr_bed}.sorted" "$cdr_bed"

log "Extracting active-HOR reads."
samtools view -@ "$threads" -b -q "$mapq" -L "$active_bed" "$bam" -o "$active_bam"
samtools index -@ "$threads" "$active_bam"
require_file "$active_bam"

log "Extracting CpG calls and inferred nucleosome footprints."
ft extract -t "$threads" --reference  --cpg "$cpg_bed12" --nuc "$nuc_bed12" "$active_bam"
bed12_blocks_to_bed6 "$cpg_bed12" "$single_cpg_bed"
bed12_blocks_to_bed6 "$nuc_bed12" "$single_nuc_bed"
require_file "$single_cpg_bed"
require_file "$single_nuc_bed"

log "Counting reference and methylated CpGs within each footprint."
bedtools nuc -fi "$reference_for_tools" -bed "$single_nuc_bed" -pattern CG |
  awk -v asm="$assembly" '
    BEGIN { OFS="\t" }
    NR>1 { print $1,$2,$3,$4,$NF,asm }
  ' | sort_bed > "$nuc_cpg_content"
require_file "$nuc_cpg_content"

bedtools intersect -sorted -a "$nuc_cpg_content" -b "$single_cpg_bed" -wao |
  awk 'BEGIN{OFS="\t"} $NF>0 && $4==$10' |
  bedtools groupby -g 1,2,3,4,5,6 -c 7 -o count > "$methylated_counts"

{
  printf 'chr\tstart\tend\tnucleosome_size_bp\tncg\tmcg\tasm\tcell\tread_id\n'
  awk -v cell="$cell" '
    BEGIN { OFS="\t" }
    FNR==NR {
      key=$1 SUBSEP $2 SUBSEP $3 SUBSEP $4
      methylated[key]=$7
      next
    }
    {
      key=$1 SUBSEP $2 SUBSEP $3 SUBSEP $4
      mcg=(key in methylated ? methylated[key] : 0)
      print $1,$2,$3,$3-$2,$5,mcg,$6,cell,$4
    }
  ' "$methylated_counts" "$nuc_cpg_content"
} > "$final_tsv"

bgzip -f "$final_tsv"
require_file "$final_gz"

if [[ "$keep_intermediates" != true ]]; then
  rm -f -- "$active_bam" "${active_bam}.bai" "$cpg_bed12" "$nuc_bed12"     "$single_cpg_bed" "$single_nuc_bed" "$nuc_cpg_content" "$methylated_counts"
fi

log "Completed successfully: $final_gz"

