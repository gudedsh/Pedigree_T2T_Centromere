#!/bin/bash
set -euo pipefail

matbam=$(realpath "$1")
patbam=$(realpath "$2")
pred_bam=$(realpath "$3")
sample="$4"
outdir=$(realpath "$5")
threads="${6:-8}"
min_len="${7:-1000}"
mg_cutoff="${8:-0}"

mkdir -p "$outdir"
cd "$outdir"

extract_score() {
    input_bam="$1"
    hap="$2"

    output_tsv="${sample}.${hap}.haplotype_scores.tsv"

    echo -e "read_id\tread_length\tmg\tNM\tnormalized_NM\trq" > "$output_tsv"

    samtools view -@ "$threads" -F 0x100 -F 0x800 "$input_bam" | \
    awk 'BEGIN{OFS="\t"}
    {
        read_id = $1
        read_len = length($10)

        mg = "NA"
        nm = "NA"
        rq = "NA"

        for (i = 12; i <= NF; i++) {
            if ($i ~ /^mg:f:/) {
                split($i, arr, ":")
                mg = arr[3]
            } else if ($i ~ /^NM:i:/) {
                split($i, arr, ":")
                nm = arr[3]
            } else if ($i ~ /^rq:f:/) {
                split($i, arr, ":")
                rq = arr[3]
            }
        }

        normalized_nm = (nm != "NA" && read_len > 0) ? nm / read_len : "NA"
        print read_id, read_len, mg, nm, normalized_nm, rq
    }' >> "$output_tsv"

    gzip -f "$output_tsv"
}

echo "[$(date)] Extracting Mat scores"
extract_score "$matbam" "Mat"

echo "[$(date)] Extracting Pat scores"
extract_score "$patbam" "Pat"

mat_file="${sample}.Mat.haplotype_scores.tsv.gz"
pat_file="${sample}.Pat.haplotype_scores.tsv.gz"

mat_tmp="Mat.score.tmp.tsv"
pat_tmp="Pat.score.tmp.tsv"
merged_tmp="merged.MatPat.score.tmp.tsv"

zcat "$mat_file" | \
awk -v min_len="$min_len" 'BEGIN{OFS="\t"} NR>1 && $2 >= min_len {print $1, $3, $5}' | \
sort -k1,1 > "$mat_tmp"

zcat "$pat_file" | \
awk -v min_len="$min_len" 'BEGIN{OFS="\t"} NR>1 && $2 >= min_len {print $1, $3, $5}' | \
sort -k1,1 > "$pat_tmp"

join -a 1 -a 2 -e "NA" -o 0,1.2,1.3,2.2,2.3 "$mat_tmp" "$pat_tmp" > "$merged_tmp"

out="phased_reads.tsv"
stats="phased_reads.stats.tsv"

echo -e "read_id\tmpphase" > "$out"

awk -v mg_cutoff="$mg_cutoff" 'BEGIN{OFS="\t"; srand(1)}
{
    read=$1
    mgM=$2
    nnmM=$3
    mgP=$4
    nnmP=$5

    phase="Unphased"

    if (mgP=="NA" && mgM!="NA") {
        phase="Mat_unique"
    } else if (mgM=="NA" && mgP!="NA") {
        phase="Pat_unique"
    } else if (mgM!="NA" && mgP!="NA") {
        del_mg = mgM - mgP
        del_nnm = nnmM - nnmP

        if (del_mg > mg_cutoff && del_nnm < 0) {
            phase="Mat_score"
        } else if (del_mg < -mg_cutoff && del_nnm > 0) {
            phase="Pat_score"
        } else {
            if (rand() < 0.5) {
                phase="Mat_random"
            } else {
                phase="Pat_random"
            }
        }
    }

    print read, phase
}' "$merged_tmp" >> "$out"

awk -v sample="$sample" 'BEGIN{OFS="\t"}
NR==1{next}
{
    total++
    split($2,a,"_")
    hap=a[1]
    type=a[2]

    if(hap=="Mat") mat++
    if(hap=="Pat") pat++
    if(type=="random") random++
    else confident++
}
END{
    print "sample","reads_N","reads_Mat","reads_Pat","reads_random","reads_confident","confident_rate"
    if(total==0){
        print sample,0,0,0,0,0,"NA"
    } else {
        print sample,total,mat+0,pat+0,random+0,confident+0,confident/total
    }
}' "$out" > "$stats"

awk 'BEGIN{FS=OFS="\t"} NR>1 && $2 ~ /^Mat_/ {print $1}' "$out" > Mat_reads_id.txt
awk 'BEGIN{FS=OFS="\t"} NR>1 && $2 ~ /^Pat_/ {print $1}' "$out" > Pat_reads_id.txt

samtools view -@ "$threads" -b -h -N Mat_reads_id.txt "$pred_bam" -o "${sample}.MatReads.bam"
samtools view -@ "$threads" -b -h -N Pat_reads_id.txt "$pred_bam" -o "${sample}.PatReads.bam"

samtools index -@ "$threads" "${sample}.MatReads.bam"
samtools index -@ "$threads" "${sample}.PatReads.bam"

rm -f "$mat_tmp" "$pat_tmp" "$merged_tmp"

echo "[$(date)] Step 3 done."
echo "Output:"
echo "  $out"
echo "  $stats"
echo "  ${sample}.MatReads.bam"
echo "  ${sample}.PatReads.bam"
