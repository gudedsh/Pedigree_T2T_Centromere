# Module 07: Call De Novo / Somatic Mutation Pipeline

This directory contains the complete computational workflow for identifying high-confidence line-specific de novo or somatic mutations using multi-platform variant calling architectures on long-read PacBio datasets. To achieve maximum baseline specificity and suppress complex structural alignment artifacts, this module integrates parallel executions of **Google DeepVariant** and **GATK4 Mutect2**, applies multi-parameter quality metrics filtering, and intersects the callsets to yield high-confidence consensus variants.

---

## Workflow Overview

The operational architecture inside this folder is distributed across three consecutive procedural stages:


### Stage 1: Variant Calling
* **`run_deepvariant.sh`**: Dispatches native multi-sharded variant calling leveraging DeepVariant's `PACBIO` neural networks to capture kinetic-aware variants.
* **`run_gatk_mutect2.sh`**: Executes traditional statistical paired somatic variant calling via GATK4 Mutect2, automatically extracting Tumor/Normal `@RG SM` tokens to screen out matching background variations.

### Stage 2: Independent Post-Filtering & Artifact Exclusions
* **`run_DV_vcf_filter.sh`**: Filters DeepVariant outputs based on strict genotype quality, allelic coverage boundaries ($10 \le \text{DP} \le 100$), and VAF limits. It strips variants overlapping low-quality Flagger zones, NucFreq anomalies, and dense Tandem Repeats (TRF).
* **`run_mutect2_vcf_filter.sh`**: Performs multi-sample GATK sorting and checks somatic indicators (`TLOD >= 6`, `NLOD >= 2`). It applies math filtering to mandate homozygous wild-type backgrounds in Normal controls while validating high-frequency somatic leakage ($> 20\%$) in Tumor tracks before running identical repeat masks exclusions.

### Stage 3: Callset Intersection & Functional Annotation
* **`intersect_dv_gatk_variants.sh`**: Intersects the clean BED coordinates from both calling workflows. Variants verified by both architectures are retained as the final consensus dataset, which is annotated across exons, promoters, centromeres, and segmental duplications (`biser`) to update `summary_final_overlap.txt`.

---

## Prerequisites & Dependencies

Ensure the following genomic toolsets are loaded in your local HPC system environment:
* **DeepVariant** (Native execution environment)
* **GATK4 Suite**
* **Bcftools**
* **Bedtools / Samtools**

---

## Usage Guide

### Run DeepVariant Callset Generation
bash run_deepvariant.sh --bam tumor.bam --ref reference.fa --outdir ./dv_raw --threads 16

### Run GATK Mutect2 Paired Calling
bash run_gatk_mutect2.sh --tumor-bam tumor.bam --normal-bam matched_pbmc.bam --ref reference.fa --outdir ./gatk_raw --threads 16

### Filter DeepVariant Outputs
bash run_DV_vcf_filter.sh --vcf ./dv_raw/sample.vcf.gz --anno-dir /path/to/annotations --outdir ./dv_filtered

### Filter GATK Mutect2 Outputs
bash run_mutect2_vcf_filter.sh --vcf ./gatk_raw/sample.vcf.gz --anno-dir /path/to/annotations --outdir ./gatk_filtered
### Final Cross-Validation Intersect
bash intersect_dv_gatk_variants.sh \
  --dv-bed ./dv_filtered/sample_clean_dv.bed \
  --gatk-bed ./gatk_filtered/sample_clean_gatk.bed \
  --anno-dir /path/to/annotations \
  --outdir ./final_consensus


