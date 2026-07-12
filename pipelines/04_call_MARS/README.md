# Module 04: Call MARS (Methylation Area Score)

This module implements **MARS (Methylation Area Score)**, a quantitative metric developed in our study to measure the magnitude of DNA hypomethylation of CDRs relative to active Higher-Order Repeat (HOR) domains.

Rather than relying solely on average methylation levels, MARS integrates both the **depth** and **genomic extent** of methylation depletion, providing a robust quantitative measurement of centromere methylation dynamics across different cell types and sequencing platforms.

The pipeline accepts genome-wide CpG methylation tracks (e.g., PacBio HiFi, Fiber-seq, or Oxford Nanopore methylation calls) together with CDR and active HOR annotations, and automatically extracts centromeric regions, performs quality control, calculates MARS scores, and generates publication-ready visualizations.

---

## Method Overview

The workflow consists of four major steps:

### 1. Extract active HOR regions

Methylation records overlapping active HOR domains and their flanking regions are efficiently extracted from bgzip-compressed methylation files using **tabix**.

### 2. Construct methylation profiles

CpG methylation values are aggregated into fixed genomic bins (default: **1 kb**) to generate continuous methylation profiles across each active HOR region and its surrounding background.

### 3. Estimate background methylation

Background methylation is estimated from regions outside the CDR after robust Winsorization, reducing the influence of local outliers and sequencing coverage fluctuations.

### 4. Calculate MARS

For each genomic bin within the active HOR,

$begin:math:display$
\\Delta\_i \= M\_\{\\mathrm\{background\}\} \- M\_i
$end:math:display$

where

- $begin:math:text$M\_i$end:math:text$ is the methylation level of the *i*-th genomic bin inside the active HOR.
- $begin:math:text$M\_\{\\mathrm\{background\}\}$end:math:text$ is the robust background methylation estimate.

The overall MARS is calculated as

$begin:math:display$
\\mathrm\{MARS\}\=\\sum\_\{i\\in \\mathrm\{HOR\}\}\\Delta\_i
$end:math:display$

The pipeline reports:

- Raw MARS
- log₂-transformed MARS
- Background methylation
- Active HOR methylation
- QC figures

---

## Workflow

```text
Compressed methylation file (.bed.gz)
                │
                ▼
        tabix region extraction
                │
                ▼
     Fixed-bin methylation profile
                │
                ▼
 Background methylation estimation
        (Robust Winsorization)
                │
                ▼
         MARS calculation
                │
                ▼
      QC plots & summary tables
```

---

## Dependencies

### External software

- Bedtools (≥2.26)
- Samtools / Tabix
- R (≥4.0)

### Required R packages

- data.table
- ggplot2
- zoo
- optparse
- scales

---

## Usage

```bash
bash call_MARS_v1.0.sh \
    -ref reference.fa \
    -cen active_hor.bed \
    -me methylation_manifest.tsv \
    -o output_directory
```

### Optional parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-bin` | 1000 | Bin size (bp) |
| `-n` | 3 | Minimum CpGs per bin |
| `-dl` | 5 | Minimum read depth |
| `-dh` | 100 | Maximum read depth |
| `-s` | 10 | Flanking distance (kb) |
| `-exp` | 0 | Active HOR expansion size (bp) |
| `-out_CDR_cutoff` | 0.10 | Winsorization cutoff outside active HOR |
| `-in_CDR_cutoff` | 0.05 | Winsorization cutoff inside active HOR |

---

## Output

For each sample, the pipeline generates:

- Processed methylation matrices
- Active HOR methylation profiles
- Background methylation estimates
- MARS summary tables
- QC figures (PDF)

---
