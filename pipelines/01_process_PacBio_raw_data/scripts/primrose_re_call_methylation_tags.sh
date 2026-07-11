#!/bin/bash

if [[ $# -lt 5 ]]; then
    echo "Usage: $0 <input_bam> <threads> <tech> <out_prefix> <test_mode_true_false>"
    exit 1
fi

BAM="$1"
THREADS="$2"
TECH="$3"
PREFIX="$4"

SUMMARY_TXT="${PREFIX}_primrose_summary.txt"
TMP_BAM="${PREFIX}_remove_meth_tags.tmp.bam"
FINAL_BAM="${PREFIX}_primrose.bam"

echo "## Check raw bam file tags ----------" > "$SUMMARY_TXT"

# 1. Dual-Verification Logic to evaluate if the origin caller is primrose
echo "## Executing tool qualification checks ----------" >> "$SUMMARY_TXT"

# Check 2.1: Audit the BAM header for the tool name registry "primrose"
HEADER_CHECK=$(samtools view -H "$BAM" | grep -i -E "^@PG.*(ID|PN):primrose" | head -n 1)

# Check 2.2: Scan the initial 1000 alignments for the specific 'fp' tag signature
TAG_STATS=$(samtools view -h "$BAM" | head -n 1000 | awk '
BEGIN {
    FS="\t"; 
    split("MM ML Mm Ml MP fp fi ri", tags, " "); 
    for (t in tags) count[tags[t]]=0
} 
!/^@/ { 
    for (i=12; i<=NF; i++) {
        split($i, a, ":"); 
        if(a[1] in count) count[a[1]]++
    }
}
END {
    for (t in count) print t ":" count[t]
}')

FP_COUNT=$(echo "$TAG_STATS" | awk -F: '$1=="fp"{print $2}')
FP_COUNT="${FP_COUNT:-0}"

# Diagnostic arbitration
IS_PRIMROSE="false"
if [[ -n "$HEADER_CHECK" ]]; then
    echo "[INFO] Confirmed primrose usage via BAM header line: $HEADER_CHECK" | tee -a "$SUMMARY_TXT"
    IS_PRIMROSE="true"
elif [[ "$FP_COUNT" -gt 0 ]]; then
    echo "[INFO] BAM header lacks metadata, but active 'fp' tags ($FP_COUNT counts) were detected." | tee -a "$SUMMARY_TXT"
    IS_PRIMROSE="true"
else
    echo "[INFO] No primrose footprints detected in header or alignment tag blocks." | tee -a "$SUMMARY_TXT"
fi

# 3. Decision Tree: Execute primrose recalculation or bypass
if [[ "$IS_PRIMROSE" == "true" ]]; then
    echo "[INFO] BAM is verified as primrose origin. Bypassing re-call phase." | tee -a "$SUMMARY_TXT"
    CURRENT_WORKING_BAM="$WORKING_INPUT_BAM"
else
    echo "[INFO] Initiating tag extraction and re-calling with primrose..." | tee -a "$SUMMARY_TXT"
    echo "## Removing outdated/non-standard methylation tags ----------" >> "$SUMMARY_TXT"
    
    # Strip previous conflicting methylation fields
    samtools view -@ "$THREADS" -h "$BAM" | \
    samtools view -@ "$THREADS" -b \
        --remove-tag MM --remove-tag ML --remove-tag MP --remove-tag Mm --remove-tag Ml \
        -o "$TMP_BAM"
    
    # Configure parameters based on technology
    PRIMROSE_OPTS="-j ${THREADS}"
    if [[ "$TECH" != "PB" ]]; then
        # Fiber-seq/DiMeLo-seq require kinetics data to allow downstream ft predict m6a modeling
        PRIMROSE_OPTS="${PRIMROSE_OPTS} --keep-kinetics"
    fi
    
    echo "Running command: primrose ${PRIMROSE_OPTS} --log-level DEBUG ${TMP_BAM} ${PREFIX}_recalled.bam" >> "$SUMMARY_TXT"
    primrose ${PRIMROSE_OPTS} --log-level DEBUG "$TMP_BAM" "${PREFIX}_recalled.bam" --log-file "${PREFIX}_primrose.log"
    
    CURRENT_WORKING_BAM="${PREFIX}_recalled.bam"
    rm -f "$TMP_BAM"
fi

########################################
## 4. Data Slimming Phase: Kinetics Evacuation
########################################
echo "## Starting final structural compaction and tag slimming ----------" >> "$SUMMARY_TXT"

if [[ "$TECH" == "PB" ]]; then
    # Classic PacBio completely strips raw kinetics fields to maximize storage efficiency and alignment speeds
    samtools view -@ "$THREADS" -h "$CURRENT_WORKING_BAM" | \
    samtools view -@ "$THREADS" -b --remove-tag fi --remove-tag ri --remove-tag fp -o "$FINAL_BAM"
    
    # Clean temporary testing files and local intermediate objects
    [[ "$WORKING_INPUT_BAM" != "$BAM" ]] && rm -f "$WORKING_INPUT_BAM"
    if [[ "$CURRENT_WORKING_BAM" != "$WORKING_INPUT_BAM" ]]; then rm -f "$CURRENT_WORKING_BAM"; fi
else
    # DML and FS retain base kinetics data for fibertools.
    # Fibertools will purge kinetics automatically afterwards in Step 1 while saving m6a arrays.
    if [[ "$CURRENT_WORKING_BAM" != "$BAM" && "$CURRENT_WORKING_BAM" != "$WORKING_INPUT_BAM" ]]; then
        mv "$CURRENT_WORKING_BAM" "$FINAL_BAM"
        [[ "$WORKING_INPUT_BAM" != "$BAM" ]] && rm -f "$WORKING_INPUT_BAM"
    else
        if [[ "$WORKING_INPUT_BAM" != "$BAM" ]]; then
            mv "$WORKING_INPUT_BAM" "$FINAL_BAM"
        else
            ln -sf "$BAM" "$FINAL_BAM"
        fi
    fi
fi

echo "[INFO] Primrose qualification phase resolved successfully." >> "$SUMMARY_TXT"
