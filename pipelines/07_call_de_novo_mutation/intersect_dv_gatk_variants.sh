#!/bin/bash

# Activate strict error trapping switches
set -euo pipefail

# Print clean usage helper instructions
usage() {
    echo "Usage: bash $0 --dv-bed <dv_output.bed> --gatk-bed <gatk_output.bed> --anno-dir <anno_dir> --outdir <output_dir>"
    echo ""
    echo "Required Arguments:"
    echo "  --dv-bed    Path to the filtered DeepVariant clean BED file."
    echo "  --gatk-bed  Path to the filtered GATK Mutect2 clean BED file."
    echo "  --anno-dir  Directory containing individual genome annotation BED/BEDPE tracks."
    echo "  --outdir    Target workspace output directory where consensus assets will be saved."
    exit 1
}

# Assign flexible fallback options variables
DV_BED=""
GATK_BED=""
ANNO_DIR=""
OUTDIR=""

# Parse positional options parameters
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dv-bed)    DV_BED="$2"; shift 2 ;;
        --gatk-bed)  GATK_BED="$2"; shift 2 ;;
        --anno-dir)  ANNO_DIR="$2"; shift 2 ;;
        --outdir)    OUTDIR="$2"; shift 2 ;;
        -h|--help)   usage ;;
        *) echo "[ERROR] Unknown option parameter argument: $1"; usage ;;
    esac
done

# Perform parameter requirement validations
if [ -z "$DV_BED" ] \vert{}\vert{} [ -z "$GATK_BED" ] || [ -z "$ANNO_DIR" ] \vert{}\vert{} [ -z "$OUTDIR" ]; then
    echo "[ERROR] Missing mandatory structural input parameters."
    usage
fi

# Resolve execution tracks to concrete absolute paths
DV_ABS=$(realpath "$DV_BED")
GATK_ABS=$(realpath "$GATK_BED")
ANNO_ABS=$(realpath "$ANNO_DIR")
mkdir -p "$OUTDIR"
OUTDIR_ABS=$(realpath "$OUTDIR")

# Deriving clean sample prefixes from the filenames
DV_BASE=$(basename "$DV_ABS")
PREFIX=$(echo "$DV_BASE" | sed -E 's/(_deepvaraints)?(_gatk)?(_somatic_rm_tandemR_Gap)?\.bed//')

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

# Verify that essential annotation files exist
for file in "$DUP" "$CT" "$GP" "$PRO" "$EXON"; do
    if [ -z "$file" ] \vert{}\vert{} [ ! -f "$file" ]; then
        echo "[ERROR] Mandatory annotation track file missing or unresolvable in $ANNO_ABS."
        exit 1
    fi
done

# Define final output path targets
FINAL_BED="${OUTDIR_ABS}/${PREFIX}_consensus_overlap_variants.bed"
SUMMARY_FILE="${OUTDIR_ABS}/summary_final_overlap.txt"

echo "========================================================="
echo " [INFO] Initializing Callset Cross-Validation Overlap"
echo " [TIME] $(date)"
echo "========================================================="
echo ">> DeepVariant BED: $DV_ABS"
echo ">> GATK Mutect2   : $GATK_ABS"
echo ">> Sample Prefix  : $PREFIX"
echo ">> Target Cohort  : $IND ($HAP)"
echo ">> Consensus Output: $FINAL_BED"

# 1. Run coordinate intersect mapping to find overlapping consensus variants
echo "[INFO] Intersecting DeepVariant and GATK callsets..."
intersectBed -a "$DV_ABS" -b "$GATK_ABS" \vert{} sortBed -i > "$FINAL_BED"

# 2. Compute absolute metric counts across partitioned genomic zones
echo "[INFO] Quantifying consensus regional overlap markers..."
ALL_COUNTS=$(cat "$FINAL_BED" | wc -l)

if [ "$ALL_COUNTS" -eq 0 ]; then
    echo "[WARN] Zero variants overlapped between both calling platforms. Logging zero arrays."
    IN_CEN=0; IN_DUP=0; IN_PRO=0; IN_EXON=0
else
    IN_CEN=$(intersectBed -a "$FINAL_BED" -b "$CT" -wao \vert{} awk '$NF>0' | cut -f 1,2,3,4,5 | sort -u | wc -l)
    IN_DUP=$(intersectBed -a "$FINAL_BED" -b "$DUP" -wao \vert{} awk '$NF>0' | cut -f 1,2,3,4,5 | sort -u | wc -l)
    IN_PRO=$(intersectBed -a "$FINAL_BED" -b "$PRO" -wao \vert{} awk '$NF>0' | cut -f 1,2,3,4,5 | sort -u | wc -l)
    IN_EXON=$(intersectBed -a "$FINAL_BED" -b "$EXON" -wao \vert{} awk '$NF>0' | cut -f 1,2,3,4,5 | sort -u | wc -l)
fi

# 3. Log consensus results to the cumulative master table
if [ ! -s "$SUMMARY_FILE" ]; then
    echo -e "Sample_ID\tConsensus_Overlap_Variants\tPromoter_Overlap\tExon_Overlap\tCentromere_Overlap\tSegDup_Overlap" > "$SUMMARY_FILE"
fi

echo -e "${PREFIX}\t${ALL_COUNTS}\t${IN_PRO}\t${IN_EXON}\t${IN_CEN}\t${IN_DUP}" >> "$SUMMARY_FILE"

echo "========================================================="
echo " [SUCCESS] Cross-validation overlap pipeline completed."
echo " [METRICS] Summary table updated inside: $SUMMARY_FILE"
echo "========================================================="


