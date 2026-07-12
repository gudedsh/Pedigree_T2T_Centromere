# Module 07: Call De Novo / Somatic Mutations

This module implements a consensus variant-calling pipeline for identifying high-confidence **de novo** and **somatic mutations** from long-read PacBio HiFi sequencing data.

To maximize specificity while minimizing false-positive calls arising from sequencing or alignment artifacts, the pipeline combines two complementary variant callers (**DeepVariant** and **GATK Mutect2**), applies stringent quality filtering, and retains only consensus variants supported by both methods.

The final variant set is further annotated with genomic features, including genes, promoters, centromeres, tandem repeats, and segmental duplications, facilitating downstream analyses of mutation distribution across genomic contexts.

---

## Method Overview

The workflow consists of three major steps.

### 1. Variant calling

Candidate variants are independently identified using two complementary approaches:

- **DeepVariant** for high-accuracy germline variant calling on PacBio HiFi reads.
- **GATK Mutect2** for paired somatic variant calling using matched normal samples.

Running both callers improves sensitivity while allowing cross-validation of candidate mutations.

### 2. Quality filtering

Variant calls from each caller are independently filtered using stringent quality criteria, including sequencing depth, genotype quality, variant allele frequency, and caller-specific confidence metrics.

To reduce false positives, variants overlapping low-confidence genomic regions (e.g., tandem repeats, low-quality assembly regions, and other predefined genomic masks) are removed.

### 3. Consensus variant identification

Filtered DeepVariant and Mutect2 callsets are intersected to retain only variants supported by both methods, producing a high-confidence consensus mutation set.

The final variants are annotated with genomic features, including:

- Genes
- Promoters
- Centromeres
- Segmental duplications (BISER)
- Other user-defined genomic annotations

---

## Workflow

```text
PacBio HiFi BAM
        │
        ├──────────────► DeepVariant
        │
        └──────────────► GATK Mutect2
                 │
                 ▼
      Independent quality filtering
                 │
                 ▼
      Consensus variant intersection
                 │
                 ▼
       Functional annotation
                 │
                 ▼
 High-confidence de novo / somatic variants
```

---

## Dependencies

### External software

- DeepVariant
- GATK4
- BCFtools
- Bedtools
- Samtools

---

## Usage

### 1. Run DeepVariant

```bash
bash run_deepvariant.sh \
    --bam tumor.bam \
    --ref reference.fa \
    --outdir dv_raw \
    --threads 16
```

### 2. Run Mutect2

```bash
bash run_gatk_mutect2.sh \
    --tumor-bam tumor.bam \
    --normal-bam matched_pbmc.bam \
    --ref reference.fa \
    --outdir gatk_raw \
    --threads 16
```

### 3. Filter DeepVariant calls

```bash
bash run_DV_vcf_filter.sh \
    --vcf dv_raw/sample.vcf.gz \
    --anno-dir annotation_directory \
    --outdir dv_filtered
```

### 4. Filter Mutect2 calls

```bash
bash run_mutect2_vcf_filter.sh \
    --vcf gatk_raw/sample.vcf.gz \
    --anno-dir annotation_directory \
    --outdir gatk_filtered
```

### 5. Generate consensus variants

```bash
bash intersect_dv_gatk_variants.sh \
    --dv-bed dv_filtered/sample_clean_dv.bed \
    --gatk-bed gatk_filtered/sample_clean_gatk.bed \
    --anno-dir annotation_directory \
    --outdir final_consensus
```

---

## Output

For each sample, the pipeline generates:

- Filtered DeepVariant callset
- Filtered Mutect2 callset
- High-confidence consensus variants
- Functional annotation tables
- Mutation summary reports

---
