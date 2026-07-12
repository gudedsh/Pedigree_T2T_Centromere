#!/bin/bash

# Activate strict error trapping switches
set -euo pipefail

# Print clean usage helper instructions
usage() {
    echo "Usage: bash $0 --vcf <input_gatk.vcf.gz> --anno-dir <anno_dir> --outdir <output_dir>"
    echo ""
    echo "Required Arguments:"
    echo "  --vcf       Absolute or relative path to the paired GATK Mutect2 VCF (.vcf.gz)."
    echo "  --anno-dir  Directory containing individual genome annotation BED/BEDPE tracks."
    echo "  --outdir    Target workspace output directory where outputs will be saved."
    exit 1
}

# Assign options variables
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

# Deriving sample prefix tokens
VCF_BASE=$(basename "$VCF_ABS")
PREFIX=$(echo "$VCF_BASE" | sed -E 's/(_aligned)?(_sorted)?(_gatk_Mutect2)?\.vcf\.gz//')

# Extract pedigree cohort descriptors (e.g., PAN010, Mat/Pat)
IND=$(echo "$PREFIX" | grep -o 'PAN[0-9]*' || echo "PAN010")
HAP=$(echo "$PREFIX" | grep -o '[MP]at' || head -n 1 || echo "Mat")

# Setup absolute file paths for annotation tracks
DUP="${ANNO_ABS}/${IND}${HAP}_biser.bedpe"
if [ ! -f "$DUP" ]; then
    DUP=$(ls ${ANNO_ABS}/*${IND}*${HAP}*_biser.bedpe 2>/dev/null \vert{} head -n 1 \vert{}\vert{} echo "${ANNO_ABS}/PAN010Mat_biser.bedpe")
fi

CT=$(ls${ANNO_ABS}/*${IND}*${HAP}*cenSat.bed 2>/dev/null | head -n 1 || echo "")
GP=$(ls${ANNO_ABS}/*${IND}*${HAP}*flagger_final.no_Hap.bed 2>/dev/null | head -n 1 || echo "")
PRO=$(ls${ANNO_ABS}/*${IND}*${HAP}*_gene_annotations_allGenes_promoter_up1kdown100bp.bed 2>/dev/null | head -n 1 || echo "")
EXON=$(ls${ANNO_ABS}/*${IND}*${HAP}*_gene_annotations_unique_exon.bed 2>/dev/null | head -n 1 || echo "")
TR=$(ls${ANNO_ABS}/*${IND}*${HAP}*mono+di.trf.minL10.slop5.bed 2>/dev/null | head -n 1 || echo "")
NFG=$(ls${ANNO_ABS}/*${IND}*${HAP}*_nucflag_conf.bed 2>/dev/null | head -n 1 || echo "")

# Verify that essential annotation files exist
for file in "$DUP" "$CT" "$GP" "$PRO" "$EXON" "$TR" "$NFG"; do
    if [ -z "$file" ] \vert{}\vert{} [ ! -f "$file" ]; then
        echo "[ERROR] Mandatory annotation track file missing or unresolvable in $ANNO_ABS."
        exit 1
    fi
done

# Define final output path targets
OUT_BEDBird="${OUTDIR_ABS}/${PREFIX}_gatk_somatic_rm_tandemR_Gap.bed"
SUMMARY_FILE="${OUTDIR_ABS}/summary_gatk.txt"
TEMP_VCF="${OUTDIR_ABS}/${PREFIX}.tem_filter.vcf.gz"

echo "========================================================="
echo " [INFO] Initializing GATK Mutect2 VCF Post-Filtering"
echo " [TIME] $(date)"
echo "========================================================="
echo ">> Input GATK VCF : $VCF_ABS"
echo ">> Sample Prefix  : $PREFIX"
echo ">> Output Clean   : $OUT_BEDBird"

# 1. Dynamically parse multi-sample headers for Tumor and Normal IDs
echo "[INFO] Resolving sample headers via bcftools..."
NORMAL_ID=$(bcftools query -l "$VCF_ABS" | grep "PBMC" | head -n 1 || echo "")
TUMOR_ID=$(bcftools query -l "$VCF_ABS" | grep -E "iPSC|NPC" | head -n 1 || echo "")

if [ -z "$NORMAL_ID" ] \vert{}\vert{} [ -z "$TUMOR_ID" ]; then
    echo "[ERROR] Failed to extract paired PBMC (Normal) and iPSC/NPC (Tumor) sample tokens from VCF header."
    exit 1
fi

echo " -> Normal Channel Sample ID: $NORMAL_ID"
echo " -> Tumor Channel Sample ID : $TUMOR_ID"

# 2. Extract and index ordered sample sub-matrices
echo "[INFO] Building structured intermediate target array..."
bcftools view -s "${NORMAL_ID},${TUMOR_ID}" "$VCF_ABS" -Oz -o "$TEMP_VCF"
bcftools index "$TEMP_VCF"

# 3. Apply mathematical VAF filters, hard thresholds, and mask exclusions
echo "[INFO] Running multi-layer somatic variant screening and bedtools intersections..."
bcftools filter "$TEMP_VCF" \
  -i "FORMAT/DP[1] >= 10 && FORMAT/AD[1:1] / FORMAT/DP[1] >= 0.1 && INFO/TLOD >= 6 && INFO/NLOD >= 2 && FORMAT/AD[0:1] / FORMAT/DP[0] < 0.05" \
  --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
  | bcftools query -f '%CHROM\t%POS\t%END\t%TYPE\t%REF\t%ALT\t[%GT\t%AD\t]\n' \
  | sed 's/,/\t/g' \
  | awk '$8+$9 >= 10 && $8+$9 <= 100 && $11+$12 >= 10 && $11+$12 <= 100 && $9/($8+$9) < 0.05 && $12/($11+$12) > 0.2 && ($7=="0|0" || $7=="0/0") {OFS="\t"; if ($4=="SNP") $3=$3+1; print $1,$2,$3,$4,$5,$6,$10,$12/($11+$12)}' \
  | sortBed -i \
  | intersectBed -a - -b "$TR" -v \
  | intersectBed -a - -b "$GP" -v \
  | intersectBed -a - -b "$NFG" -v > "$OUT_BEDBird"

# Clean up intermediate tracking matrix files
rm -f "$TEMP_VCF" "${TEMP_VCF}.csi"

# 4. Compute structural regional overlap counts
echo "[INFO] Quantifying structural regional overlap markers..."
ALL_COUNTS=$(cat "$OUT_BEDBird" | wc -l)

if [ "$ALL_COUNTS" -eq 0 ]; then
    echo "[WARN] Zero variants passed strict artifact cleanup filtering."
    IN_CEN=0; IN_DUP=0; IN_PRO=0; IN_EXON=0
else
    IN_CEN=$(intersectBed -a "$OUT_BEDBird" -b "$CT" -wao \vert{} awk '$NF>0' | cut -f 1,2,3,4,5 | sort -u | wc -l)
    IN_DUP=$(intersectBed -a "$OUT_BEDBird" -b "$DUP" -wao \vert{} awk '$NF>0' | cut -f 1,2,3,4,5 | sort -u | wc -l)
    IN_PRO=$(intersectBed -a "$OUT_BEDBird" -b "$PRO" -wao \vert{} awk '$NF>0' | cut -f 1,2,3,4,5 | sort -u | wc -l)
    IN_EXON=$(intersectBed -a "$OUT_BEDBird" -b "$EXON" -wao \vert{} awk '$NF>0' | cut -f 1,2,3,4,5 | sort -u | wc -l)
fi

# 5. Log results to the cumulative master table
if [ ! -s "$SUMMARY_FILE" ]; then
    echo -e "Sample_ID\tTotal_Somatic_Variants\tPromoter_Overlap\tExon_Overlap\tCentromere_Overlap\tSegDup_Overlap" > "$SUMMARY_FILE"
fi

echo -e "${PREFIX}\t${ALL_COUNTS}\t${IN_PRO}\t${IN_EXON}\t${IN_CEN}\t${IN_DUP}" >> "$SUMMARY_FILE"

echo "========================================================="
echo " [SUCCESS] GATK Somatic filtering workflow completed."
echo " [METRICS] Summary table updated inside: $SUMMARY_FILE"
echo "========================================================="

