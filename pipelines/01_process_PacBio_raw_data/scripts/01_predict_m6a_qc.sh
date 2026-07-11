#!/bin/bash
set -euo pipefail

bam=$(realpath "$1")
threads=${2:-8}
prefix="$3"
outdir=$(realpath "$4")

mkdir -p "$outdir"
cd "$outdir"

pred_bam="${prefix}.predict_m6a.bam"
per_read="${prefix}.dimelo_m6a_per_read.tsv"
summary="${prefix}.dimelo_m6a_qc.summary.txt"
len_dist="${prefix}.read_length_distribution.tsv"
log="${prefix}.dimelo_m6a_qc.log"

echo "Input BAM: $bam" | tee "$log"
echo "Threads: $threads" | tee -a "$log"

echo "[$(date)] Running ft predict-m6a..." | tee -a "$log"

sk() {
ft --version | tee -a "$log"

ft predict-m6a \
  -t "$threads" \
  "$bam" \
  "$pred_bam" 2>&1 | tee -a "$log"

samtools index -@ "$threads" "$pred_bam"
}
echo -e "read_id\tread_len\tn_m6a\tm6a_per_kb\thas_m6a" > "$per_read"

samtools view -@ "$threads" "$pred_bam" | \
awk 'BEGIN{OFS="\t"}
{
    read_id=$1
    read_len=length($10)
    n_m6a=0

    for(i=12;i<=NF;i++){
        if($i ~ /^MM:Z:/){
            mm=$i
            sub(/^MM:Z:/,"",mm)
            n=split(mm, groups, ";")
            for(g=1; g<=n; g++){
                if(groups[g] ~ /^A[+-]a[.?]?,/){
                    split(groups[g], arr, ",")
                    n_m6a += length(arr)-1
                }
            }
        }
    }

    if(read_len > 0){
          has_m6a = (n_m6a > 0) ? 1 : 0 ; print read_id, read_len, n_m6a, n_m6a/read_len*1000, has_m6a 
    }
}' >> "$per_read"

awk 'BEGIN{FS=OFS="\t"}
NR>1{
    n++
    len[n]=$2
    total_len += $2
    total_m6a += $3
    total_density += $4
    m6a_reads += $5
}
END{
    if(n==0){print "reads",0; exit}
    print "reads", n
    print "total_bases", total_len
    print "mean_read_len", total_len/n
    print "mean_m6a_per_read", total_m6a/n
    print "total_m6a_calls", total_m6a
    print "m6a_positive_reads", m6a_reads
    print "m6a_positive_read_fraction", m6a_reads/n
    print "mean_m6a_per_kb", total_density/n
}' "$per_read" > "$summary"

tail -n +2 "$per_read" | cut -f2 | sort -n | \
awk '
{
    a[NR]=$1
    total += $1
}
END{
    if(NR==0) exit
    p10=int(NR*0.10); if(p10<1)p10=1
    p25=int(NR*0.25); if(p25<1)p25=1
    p50=int(NR*0.50); if(p50<1)p50=1
    p75=int(NR*0.75); if(p75<1)p75=1
    p90=int(NR*0.90); if(p90<1)p90=1

    print "read_len_p10", a[p10]
    print "read_len_p25", a[p25]
    print "read_len_median", a[p50]
    print "read_len_p75", a[p75]
    print "read_len_p90", a[p90]

    half=total/2
    cumsum=0
    for(i=NR;i>=1;i--){
        cumsum += a[i]
        if(cumsum >= half){
            print "read_len_N50", a[i]
            break
        }
    }
}' >> "$summary"

echo -e "\n## m6A per read distribution" >> "$summary"

tail -n +2 "$per_read" | cut -f3 | sort -n | \
awk '
{
    a[NR]=$1
}
END{
    if(NR==0) exit
    p10=int(NR*0.10); if(p10<1)p10=1
    p25=int(NR*0.25); if(p25<1)p25=1
    p50=int(NR*0.50); if(p50<1)p50=1
    p75=int(NR*0.75); if(p75<1)p75=1
    p90=int(NR*0.90); if(p90<1)p90=1

    print "m6a_per_read_p10", a[p10]
    print "m6a_per_read_p25", a[p25]
    print "m6a_per_read_median", a[p50]
    print "m6a_per_read_p75", a[p75]
    print "m6a_per_read_p90", a[p90]
}' >> "$summary"

echo -e "length_bin\treads" > "$len_dist"

tail -n +2 "$per_read" | \
awk 'BEGIN{OFS="\t"}
{
    len=$2
    if(len < 1000) bin="<1kb"
    else if(len < 2000) bin="1-2kb"
    else if(len < 5000) bin="2-5kb"
    else if(len < 10000) bin="5-10kb"
    else if(len < 15000) bin="10-15kb"
    else if(len < 20000) bin="15-20kb"
    else if(len < 30000) bin="20-30kb"
    else if(len < 50000) bin="30-50kb"
    else bin=">50kb"
    count[bin]++
}
END{
    bins[1]="<1kb"; bins[2]="1-2kb"; bins[3]="2-5kb"
    bins[4]="5-10kb"; bins[5]="10-15kb"; bins[6]="15-20kb"
    bins[7]="20-30kb"; bins[8]="30-50kb"; bins[9]=">50kb"
    for(i=1;i<=9;i++) print bins[i], count[bins[i]]+0
}' >> "$len_dist"

echo "[$(date)] Step 1 done." | tee -a "$log"

