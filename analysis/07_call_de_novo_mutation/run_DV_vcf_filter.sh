#!/bin/bash

# Activate strict error trapping switches
set -euo pipefail

# Print clean usage helper instructions
usage() {
    echo "Usage: bash $0 --vcf <input.vcf.gz> --anno-dir <anno_dir> --outdir <output_dir>"
    echo ""
    echo "Required Arguments:"
    echo "  --vcf       Absolute or relative path to the input raw DeepVariant VCF file (.vcf.gz)."
    echo "  --anno-dir  Directory containing individual genome annotation BED/BEDPE tracks."
    echo "  --outdir    Target workspace output directory where filtered tracks and summary will be saved."
    exit 1
}

# Assign flexible fallback options variables
VCF_INPUT=""
ANNO_DIR=""
OUTDIR=""

# Parse positional options parameters
while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)       VCF_INPUT="$2"; shift 2 ;;
        --anno-dir)  ANNO_DIR="$2"; shift 2 ;;
        --outdir)    OUTDIR="$2"; shift 2 ;;
        -h|--help)   usage ;;
        *) echo "[ERROR] Unknown option parameter argument: $1"; usage ;;
    esac
done

# Perform parameter requirement validations
if [ -z "$VCF_INPUT" ] || [ -z "$ANNO_DIR" ] \vert{}\vert{} [ -z "$OUTDIR" ]; then
    echo "[ERROR] Missing mandatory structural input parameters."
    usage
fi

# Resolve execution tracks to concrete absolute paths
VCF_ABS=$(realpath "$VCF_INPUT")
ANNO_ABS=$(realpath "$ANNO_DIR")
mkdir -p "$OUTDIR"
OUTDIR_ABS=$(realpath "$OUTDIR")

# Deriving sample descriptors for robust filename tracking
VCF_BASE=$(basename "$VCF_ABS")
# Automatically strips common suffix variations to isolate clean prefix
PREFIX=$(echo "$VCF_BASE" | sed -E 's/(_aligned)?(_sorted)?(_deepvaraints_callGV|_callGV)?\.vcf\.gz//')

# Extract pedigree cohort descriptors (e.g., PAN010, Mat/Pat)
IND=$(echo "$PREFIX" | grep -o 'PAN[0-9]*' || echo "PAN010")
HAP=$(echo "$PREFIX" | grep -o '[MP]at' || echo "Mat")

# Setup absolute file paths for annotation tracks
DUP="${ANNO_ABS}/${IND}${HAP}_biser.bedpe"
# Fallback logic to catch case variants in filenames (e.g., PAN010Mat vs PAN010_Mat)
if [ ! -f "$DUP" ]; then
    DUP=$(ls ${ANNO_ABS}/*${IND}*${HAP}*_biser.bedpe 2>/dev/null \vert{} head -n 1 \vert{}\vert{} echo "${ANNO_ABS}/PAN010Mat_biser.bedpe")
fi

CT=$(ls${ANNO_ABS}/*${IND}*${HAP}*cenSat.bed 2>/dev/null | head -n 1 || echo "")
GP=$(ls${ANNO_ABS}/*${IND}*${HAP}*flagger_final.no_Hap.bed 2>/dev/null | head -n 1 || echo "")
PRO=$(ls${ANNO_ABS}/*${IND}*${HAP}*_gene_annotations_allGenes_promoter_up1kdown100bp.bed 2>/dev/null | head -n 1 || echo "")
EXON=$(ls${ANNO_ABS}/*${IND}*${HAP}*_gene_annotations_unique_exon.bed 2>/dev/null | head -n 1 || echo "")
TR=$(ls${ANNO_ABS}/*${IND}*${HAP}*mono+di.trf.minL10.slop5.bed 2>/dev/null | head -n 1 || echo "")
NFG=$(ls${ANNO_ABS}/*${IND}*${HAP}*_nucflag_conf.bed 2>/dev/null | head -n 1 || echo "")

# Verify that essential annotation files exist before running downstream intersects
for file in "$DUP" "$CT" "$GP" "$PRO" "$EXON" "$TR" "$NFG"; do
    if [ -z "$file" ] \vert{}\vert{} [ ! -f "$file" ]; then
        echo "[ERROR] Mandatory annotation track file missing or unresolvable in $ANNO_ABS. Please check filename templates."
        exit 1
    fi
done

# Define final output path targets
OUT_BED="${OUTDIR_ABS}/${PREFIX}_deepvaraints_somatic_rm_tandemR_Gap.bed"
SUMMARY_FILE="${OUTDIR_ABS}/summary.txt"

echo "========================================================="
echo " [INFO] Initializing VCF Post-Filtering and Annotation"
echo " [TIME] $(date)"
echo "========================================================="
echo ">> Input VCF File   : $VCF_ABS"
echo ">> Sample Prefix    : $PREFIX"
echo ">> Individual Token : $IND ($HAP)"
echo ">> Output Clean BED : $OUT_BED"
echo ">> Target Master Log: $SUMMARY_FILE"

# 1. Execute Bcftools multi-metric filtering and Bedtools sequential clean-up exclusions
echo "[INFO] Commencing bcftools filtering and raw artifact exclusions..."
bcftools view -f PASS -V other \
  -i '(FORMAT/GT = "1/1" || FORMAT/GT = "0/1") && FORMAT/DP >= 10 && FORMAT/DP <= 100 && FORMAT/VAF > 0.1' \
  --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
  "$VCF_ABS" \
  | bcftools query -f '%CHROM\t%POS\t%END\t%TYPE\t%REF\t%ALT\t[%GT\t%VAF\t]\n' \
  | awk 'BEGIN{OFS="\t"} {if ($4=="SNP") $3=$3+1; print $1,$2,$3,$4,$5,$6,$7,$8}' \
  | intersectBed -a - -b "$TR" -v \
  | intersectBed -a - -b "$GP" -v \
  | intersectBed -a - -b "$NFG" -v > "$OUT_BED"

# 2. Compute absolute metric counts across partitioned genomic zones
echo "[INFO] Quantifying structural regional overlap markers..."
ALL_COUNTS=$(cat "$OUT_BED" | wc -l)

if [ "$ALL_COUNTS" -eq 0 ]; then
    echo "[WARN] Zero variants passed strict artifact cleanup filtering. Logging zero distribution metrics."
    IN_CEN=0; IN_DUP=0; IN_PRO=0; IN_EXON=0
else
    IN_CEN=$(intersectBed -a "$OUT_BED" -b "$CT" -wao \vert{} awk '$NF>0' | cut -f 1,2,3,4,5 | sort -u | wc -l)
    IN_DUP=$(intersectBed -a "$OUT_BED" -b "$DUP" -wao \vert{} awk '$NF>0' | cut -f 1,2,3,4,5 | sort -u | wc -l)
    IN_PRO=$(intersectBed -a "$OUT_BED" -b "$PRO" -wao \vert{} awk '$NF>0' | cut -f 1,2,3,4,5 | sort -u | wc -l)
    IN_EXON=$(intersectBed -a "$OUT_BED" -b "$EXON" -wao \vert{} awk '$NF>0' | cut -f 1,2,3,4,5 | sort -u | wc -l)
fi

# 3. Log results to the cumulative master table
# Creates headers if the matrix file is blank or does not exist yet
if [ ! -s "$SUMMARY_FILE" ]; then
    echo -e "Sample_ID\tTotal_Clean_Variants\tPromoter_Overlap\tExon_Overlap\tCentromere_Overlap\tSegDup_Overlap" > "$SUMMARY_FILE"
fi

echo -e "${PREFIX}\t${ALL_COUNTS}\t${IN_PRO}\t${IN_EXON}\t${IN_CEN}\t${IN_DUP}" >> "$SUMMARY_FILE"

echo "========================================================="
echo " [SUCCESS] Filtering workflow successfully completed."
echo " [METRICS] Summary table updated inside: $SUMMARY_FILE"
echo "========================================================="
