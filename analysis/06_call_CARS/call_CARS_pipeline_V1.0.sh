#!/bin/bash
#SBATCH --job-name=call_CARS_MARS_from_CENPA_DiMeLo
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
  bash call_CARS.sh \\
    -dml bam_files.tsv \\
    -cen active_hor_cdr.bed \\
    -ref reference.fa \\
    -o output_dir \\
    -bin 5000 \\
    -t 8 \\
    -s 10 \\
    -exp 0 \\
    -out_CDR_cutoff 0.10 \\
    -in_CDR_cutoff 0.05

Required:
  -dml     two-column TSV: sample_id realpath_aligned_bam_with_fiber_tags
  -cen    active HOR / CDR BED: chr active_start active_end cdr_start cdr_end
  -ref    reference fasta with .fai index

Optional:
  -o      output directory [default: ./output_results]
  -bin    bin size [default: 5000]
  -t      threads [default: SLURM_CPUS_PER_TASK or 8]
  -s      smoothing window size [default: 10]
  -exp    symmetric expansion size applied to both CDR boundaries [default: 0]
  -out_CDR_cutoff  winsorization fraction for background outside CDR [default: 0.10]
  -in_CDR_cutoff   winsorization fraction inside CDR [default: 0.05]
EOF
}



dml_file=""
activehor=""
ref_fa=""
outdir="output_results"
winsize=5000
threads=8
smooth_window=10
expansion_size=0
out_CDR_cutoff=0.10
in_CDR_cutoff=0.05


while [[ $# -gt 0 ]]; do
  case "$1" in
    -dml) dml_file="$2"; shift 2 ;;
    -cen) activehor="$2"; shift 2 ;;
    -ref) ref_fa="$2"; shift 2 ;;
    -o) outdir="$2"; shift 2 ;;
    -bin) winsize="$2"; shift 2 ;;
    -t) threads="$2"; shift 2 ;;
    -s) smooth_window="$2"; shift 2 ;;
    -exp) expansion_size="$2"; shift 2 ;;
    -out_CDR_cutoff) out_CDR_cutoff="$2"; shift 2 ;;
    -in_CDR_cutoff) in_CDR_cutoff="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done


step0_check_parameters() {

	source activate snakemake

  [[ -z "$dml_file" || -z "$activehor" || -z "$ref_fa" ]] && { usage; exit 1; }

  [[ -f "$dml_file" ]] || { echo "ERROR: bam file list not found: $dml_file"; exit 1; }
  [[ -f "$activehor" ]] || { echo "ERROR: active HOR BED not found: $activehor"; exit 1; }
  [[ -f "$ref_fa" ]] || { echo "ERROR: reference fasta not found: $ref_fa"; exit 1; }
  [[ -f "${ref_fa}.fai" ]] || { echo "ERROR: reference index not found: ${ref_fa}.fai"; exit 1; }

  command -v samtools >/dev/null || { echo "ERROR: samtools not found"; exit 1; }
  command -v bedtools >/dev/null || { echo "ERROR: bedtools not found"; exit 1; }
  command -v ft >/dev/null || { echo "ERROR: fibertools ft not found"; exit 1; }

  awk '/^#/ || NF == 0 {next} NF < 2 {print "ERROR: fs file must have at least 2 columns at line " NR;exit 1}' "$dml_file"

  while read -r sample bam; do
    [[ -z "${sample:-}" ]] && continue
    [[ "$sample" =~ ^# ]] && continue

    [[ -f "$bam" ]] || { echo "ERROR: BAM not found: $bam"; exit 1; }
    [[ -f "${bam}.bai" || -f "${bam%.bam}.bai" || -f "${bam}.csi" ]] || {
      echo "ERROR: BAM index not found for $bam"
      exit 1
    }
  done < "$dml_file"

  mkdir -p "$outdir"
  outdir=$(realpath "$outdir")

  workdir="${outdir}/intermediate_files"
  cache_bam_dir="${outdir}/cached_activeHOR_bam"
  mkdir -p "$workdir" "$cache_bam_dir" "${workdir}/sort_tmp"

  log="${outdir}/call_CARS.log"

  if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
    mem=$(( SLURM_MEM_PER_NODE / 1024 ))
  else
    mem=80
  fi

  SORT_MEM=$(( mem / (2 * threads) ))
  (( SORT_MEM < 1 )) && SORT_MEM=1

  SORTCMD="sort -k1,1 -k2,2n --parallel=${threads} -S ${SORT_MEM}G -T ${workdir}/sort_tmp"

  {
    echo "=========== Call FAS ============"
    date
    echo "BAM list: $dml_file"
    echo "Active HOR BED: $activehor"
    echo "Reference: $ref_fa"
    echo "Output directory: $outdir"
    echo "Bin size: $winsize"
    echo "Threads: $threads"
    echo "Memory: ${mem}G"
    echo "Sort memory per thread: ${SORT_MEM}G"
  } > "$log"
}


step1_prepare_genome_bins() {
  echo -e "\n[$(date)] [step1] Preparing active HOR bins and reference counts..." | tee -a "$log"

  active_bed="${workdir}/step1_active_HOR_regions.bed"
  binned_active="${workdir}/step1_active_HOR_bin${winsize}.bed"
  bin_stats="${workdir}/step1_active_HOR_bin${winsize}_ref_A_T_CpG_count.bed"

  awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' "$activehor" | $SORTCMD > "$active_bed"

  bedtools makewindows -b "$active_bed" -w "$winsize" | $SORTCMD > "$binned_active"

  bedtools getfasta -fi "$ref_fa" -bed "$binned_active" -tab |
  awk 'BEGIN{OFS="\t"}
  {
    seq=toupper($2)
    Acount=gsub(/A/,"A",seq)
    Tcount=gsub(/T/,"T",seq)
    CpGcount=gsub(/CG/,"CG",seq)
    split($1,a,":")
    split(a[2],b,"-")
    print a[1],b[1],b[2],Acount,Tcount,Acount+Tcount,CpGcount
  }' | $SORTCMD > "$bin_stats"

  awk '
  BEGIN{OFS="\t"}
  {
    n++; A+=$4; T+=$5; AT+=$6; CpG+=$7
  }
  END{
    if(n>0){
      printf("  Number of bins: %d\n", n)
      printf("  Mean A per bin: %.2f\n", A/n)
      printf("  Mean T per bin: %.2f\n", T/n)
      printf("  Mean A+T per bin: %.2f\n", AT/n)
      printf("  Mean CpG per bin: %.2f\n", CpG/n)
    }
  }' "$bin_stats" | tee -a "$log"
}



bed12_blocks_to_bed() {
  inbed="$1"
  outbed="$2"

  awk '
  BEGIN{OFS="\t"}
  {
    chr=$1
    start=$2
    read=$4
    strand=$6
    n=$10

    split($11,size,",")
    split($12,offset,",")

    for(i=1; i<=n; i++){
      s=start+offset[i]
      e=s+size[i]
      if(e<=s){
        e=s+1
      }
      print chr,s,e,read,".",strand
    }
  }' "$inbed" |
  $SORTCMD > "$outbed"
}

count_single_signal_per_read_bin() {
  read_bin_count="$1"
  signal_single_sorted="$2"
  outfile="$3"
  signal_name="$4"

  echo "  counting ${signal_name} per read-bin..." | tee -a "$log"

  bedtools intersect -sorted -a "$read_bin_count" -b "$signal_single_sorted" -wao |
  awk '
  BEGIN{OFS="\t"}
  {
    key=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8
    if($NF > 0 && $4 == $12) count[key]++
    else if(!(key in count)) count[key]=0
  }
  END{
    for(k in count) print k,count[k]
  }' | $SORTCMD > "$outfile"
}


count_patch_bp_per_read_bin() {
  read_bin_count="$1"
  patch_bed_sorted="$2"
  outfile="$3"
  signal_name="$4"

  echo "  counting ${signal_name} overlap bp per read-bin..." | tee -a "$log"

  bedtools intersect -sorted -a "$read_bin_count" -b "$patch_bed_sorted" -wo |
  awk '
  BEGIN{OFS="\t"}
  {
    key=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8
    bp[key]+=$NF
  }
  END{
    for(k in bp) print k,bp[k]
  }' | $SORTCMD > "$outfile"

  awk '
  BEGIN{OFS="\t"}
  FNR==NR{
    seen[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8]=$9
    next
  }
  {
    key=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8
    val=(key in seen ? seen[key] : 0)
    print $1,$2,$3,$4,$5,$6,$7,$8,val
  }' "$outfile" "$read_bin_count" | $SORTCMD > "${outfile}.tmp"

  mv "${outfile}.tmp" "$outfile"
}


step2_extract_activeHOR_bam() {
  sample="$1"
  bam="$2"

  bam=$(realpath "$bam")
  active_bam="${cache_bam_dir}/${sample}.activeHOR.bam"

  echo -e "\n[$(date)] [step2] Extract active HOR BAM: ${sample}" | tee -a "$log"

  if [[ -s "$active_bam" && \
        ( -f "${active_bam}.bai" || -f "${active_bam%.bam}.bai" || -f "${active_bam}.csi" ) ]]; then
    echo "  cached active HOR BAM exists, skipping:" | tee -a "$log"
    echo "    $active_bam" | tee -a "$log"
  else
    samtools view -@ "$threads" -b -L "$active_bed" "$bam" > "$active_bam"
    samtools index -@ "$threads" "$active_bam"
    echo "  active HOR BAM written: $active_bam" | tee -a "$log"
  fi
}


step3_process_one_sample() {
  sample="$1"
  active_bam="$2"

  echo -e "\n[$(date)] [step3] Process sample: ${sample}" | tee -a "$log"

  prefix="${workdir}/${sample}_Bin$((${winsize}/1000))kb"

  read_pos="${prefix}.step3_read_pos.bed"
  read_bin="${prefix}.step4_read_bin_fullcover.bed"
  read_bin_count="${prefix}.step5_read_bin_ref_count.bed"

  m6a_bed="${prefix}.step6_ft_m6a.bed"
  mcg_bed="${prefix}.step6_ft_mcg.bed"
  msp_bed="${prefix}.step6_ft_msp.bed"
  nuc_bed="${prefix}.step6_ft_nuc.bed"

  single_m6a_sorted="${prefix}.step7_single_m6a.sorted.bed"
  single_mcg_sorted="${prefix}.step7_single_mcg.sorted.bed"
  msp_sorted="${prefix}.step7_msp.sorted.bed"
  nuc_sorted="${prefix}.step7_nuc.sorted.bed"

  read_bin_m6a="${prefix}.step8_read_bin_m6a.bed"
  read_bin_mcg="${prefix}.step8_read_bin_mcg.bed"
  read_bin_msp="${prefix}.step8_read_bin_msp_bp.bed"
  read_bin_nuc="${prefix}.step8_read_bin_nuc_bp.bed"

  read_bin_all="${prefix}.step9_read_bin_all_features.tsv"
  bin_summary="${prefix}.step10_FAS_bin_summary.tsv"

  bedtools bamtobed -i "$active_bam" | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' |  $SORTCMD > "$read_pos"
  bedtools intersect -sorted -a "$read_pos" -b "$binned_active" -wao | awk -v bin="$winsize" 'BEGIN{OFS="\t"} $NF == bin {print $5,$6,$7,$4}' | $SORTCMD > "$read_bin"
  bedtools intersect -sorted -a "$read_bin" -b "$bin_stats" -wao | awk 'BEGIN{OFS="\t"} $NF > 0 {print $1,$2,$3,$4,$8,$9,$10,$11}' | $SORTCMD > "$read_bin_count"

  ft extract -t "$threads" --reference \
    --m6a "$m6a_bed" \
    --cpg "$mcg_bed" \
    --msp "$msp_bed" \
    --nuc "$nuc_bed" \
    "$active_bam"

  bed12_blocks_to_bed "$m6a_bed" "$single_m6a_sorted"
  bed12_blocks_to_bed "$mcg_bed" "$single_mcg_sorted"
  bed12_blocks_to_bed "$msp_bed" "$msp_sorted"
  bed12_blocks_to_bed "$nuc_bed" "$nuc_sorted"


  count_single_signal_per_read_bin "$read_bin_count" "$single_m6a_sorted" "$read_bin_m6a" "m6A"
  count_single_signal_per_read_bin "$read_bin_count" "$single_mcg_sorted" "$read_bin_mcg" "mCG"
  count_patch_bp_per_read_bin "$read_bin_count" "$msp_sorted" "$read_bin_msp" "MSP"
  count_patch_bp_per_read_bin "$read_bin_count" "$nuc_sorted" "$read_bin_nuc" "nucleosome"

  awk '
  BEGIN{OFS="\t"}
  FNR==NR{
    key=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8
    mcg[key]=$9
    next
  }
  {
    key=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,(key in mcg ? mcg[key] : 0)
  }' "$read_bin_mcg" "$read_bin_m6a" |
  awk '
  BEGIN{OFS="\t"}
  FNR==NR{
    key=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8
    msp[key]=$9
    next
  }
  {
    key=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8
    print $0,(key in msp ? msp[key] : 0)
  }' "$read_bin_msp" - |
  awk '
  BEGIN{OFS="\t"}
  FNR==NR{
    key=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8
    nuc[key]=$9
    next
  }
  {
    key=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8
    print $0,(key in nuc ? nuc[key] : 0)
  }' "$read_bin_nuc" - | awk -v sid="$sample" '
BEGIN{
  OFS="\t"
  print "chr","start","end","sample","read_id",
        "ref_A","ref_T","ref_AT","ref_CpG","m6A_count",
        "mCG_count_raw","mCG_count","mCG_qc_flag","MSP_bp","nuc_bp"
}
{
  ref_CpG = $8
  mCG_raw = $10
  mCG_corr = mCG_raw
  flag = "PASS"

  if(mCG_raw > ref_CpG + 2){
    flag = "FAIL_mCG_excess"
    mCG_corr = ref_CpG
  } else if(mCG_raw > ref_CpG){
    flag = "EDGE_mCG_plus2"
    mCG_corr = ref_CpG
  }

  print $1,$2,$3,sid,$4,
        $5,$6,$7,$8,
        $9,
        mCG_raw,mCG_corr,flag,
        $11,$12
}'   > "$read_bin_all"




awk -v sid="$sample" -v bin="$winsize" '
BEGIN{OFS="\t"}
NR==1 {next}
{
  key=$1"\t"$2"\t"$3

  ref_A[key]=$6
  ref_T[key]=$7
  ref_AT[key]=$8
  ref_CpG[key]=$9

  read_n[key]++

  A_denom[key]+=$6
  AT_denom[key]+=$8
  CpG_denom[key]+=$9
  bin_bp_denom[key]+=bin

  m6a_n[key]+=$10
  mcg_raw_n[key]+=$11
  mcg_n[key]+=$12
  msp_bp[key]+=$14
  nuc_bp[key]+=$15

  if($13 == "PASS") pass_n[key]++
  else if($13 == "EDGE_mCG_plus2") edge_n[key]++
  else if($13 == "FAIL_mCG_excess") fail_n[key]++
}
END{
  print "chr","start","end","sample",
        "ref_A","ref_T","ref_AT","ref_CpG","read_depth",
        "A_denom","AT_denom","CpG_denom","bin_bp_denom",
        "m6A_count","mCG_count_raw","mCG_count", "MSP_bp","nuc_bp",
        "mCG_PASS_n","mCG_EDGE_plus2_n","mCG_FAIL_excess_n",
        "m6A_per_A","m6A_per_AT","mCG_per_CpG",
        "MSP_fraction","nuc_fraction"

  for(k in read_n){
    split(k,a,"\t")

    m6a_A  = (A_denom[k]      > 0 ? m6a_n[k]/A_denom[k]       : "NA")
    m6a_AT = (AT_denom[k]     > 0 ? m6a_n[k]/AT_denom[k]      : "NA")
    mcg    = (CpG_denom[k]    > 0 ? mcg_n[k]/CpG_denom[k]     : "NA")
    msp_f  = (bin_bp_denom[k] > 0 ? msp_bp[k]/bin_bp_denom[k] : "NA")
    nuc_f  = (bin_bp_denom[k] > 0 ? nuc_bp[k]/bin_bp_denom[k] : "NA")

    printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t",
           a[1],a[2],a[3],sid,
           ref_A[k],ref_T[k],ref_AT[k],ref_CpG[k],
           read_n[k],
           A_denom[k],AT_denom[k],CpG_denom[k],bin_bp_denom[k],
           m6a_n[k],mcg_raw_n[k],mcg_n[k],
           msp_bp[k],nuc_bp[k],
           pass_n[k]+0,edge_n[k]+0,fail_n[k]+0)

    if(m6a_A=="NA") printf("NA\t"); else printf("%.6f\t",m6a_A)
    if(m6a_AT=="NA") printf("NA\t"); else printf("%.6f\t",m6a_AT)
    if(mcg=="NA") printf("NA\t"); else printf("%.6f\t",mcg)
    if(msp_f=="NA") printf("NA\t"); else printf("%.6f\t",msp_f)
    if(nuc_f=="NA") printf("NA\n"); else printf("%.6f\n",nuc_f)
  }
}' "$read_bin_all" > "$bin_summary"


  echo "  read-bin output: $read_bin_all" | tee -a "$log"
  echo "  bin-level summary: $bin_summary" | tee -a "$log"
}


step4_combine_results () {
  combined="${outdir}/CARS_MARS_bin_summary.tsv"

  first=1
  for f in "$workdir"/*.step10_FAS_bin_summary.tsv; do
    [[ -f "$f" ]] || continue
    if [[ "$first" -eq 1 ]]; then
      cat "$f" > "$combined"
      first=0
    else
      tail -n +2 "$f" >> "$combined"
    fi
  done

  echo -e "\nCombined CARS MARS bin summary written to: $combined" | tee -a "$log"
}



step5_call_CARS_MARS_r () {
 source activate r_env

  echo -e "\n[$(date)] [step5] Calling CARS and MARS in R..." | tee -a "$log"

  combined="${outdir}/CARS_MARS_bin_summary.tsv"
  cars_mars_out="${outdir}/call_CARS_results_from_CENPA_DiMeLo.tsv"
  cars_mars_plotdir="${outdir}/call_CARS_QC_plots"

  rcode="$(realpath ./r_code/call_cars_mars.r)"

  [[ -f "$rcode" ]] || {
    echo "ERROR: R code not found: $rcode" | tee -a "$log"
    exit 1
  }

  [[ -s "$combined" ]] || {
    echo "ERROR: combined FAS bin summary not found or empty: $combined" | tee -a "$log"
    exit 1
  }


  Rscript "$rcode" \
    --input "$combined" \
    --cen "$activehor" \
    --out "$cars_mars_out" \
    --plotdir "$cars_mars_plotdir" \
    --smooth_k "$smooth_window" \
    --bin_size "$winsize" \
    --cars_col "m6A_per_AT" \
    --out_CDR_cutoff "$out_CDR_cutoff" \
    --in_CDR_cutoff "$in_CDR_cutoff" \
    --cdr_expansion "$expansion_size" \
    >> "$log" 2>&1

  echo "CARS/MARS results written to: $cars_mars_out" | tee -a "$log"
}



step0_check_parameters

step1_prepare_genome_bins

while read -r sample bam; do
  [[ -z "${sample:-}" ]] && continue
  [[ "$sample" =~ ^# ]] && continue
  step2_extract_activeHOR_bam "$sample" "$bam"
done < "$dml_file"

while read -r sample bam; do
  [[ -z "${sample:-}" ]] && continue
  [[ "$sample" =~ ^# ]] && continue
  active_bam="${cache_bam_dir}/${sample}.activeHOR.bam"
  step3_process_one_sample "$sample" "$active_bam"
done < "$dml_file"

step4_combine_results

step5_call_CARS_MARS_r

{
  echo "All done."
  date
} | tee -a "$log"






