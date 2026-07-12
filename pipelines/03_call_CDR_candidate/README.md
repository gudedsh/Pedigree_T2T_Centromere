# Module 03: Call CDR (Centromere Dip Region) Candidates

This module identifies **Centromere Dip Regions (CDRs)**, localized hypomethylated domains within active Higher-Order Repeat (HOR) arrays.

The pipeline accepts bgzip-compressed CpG methylation tracks together with active HOR annotations. It automatically extracts methylation records from active HOR regions, constructs binned methylation profiles, identifies candidate CDRs, performs quality control, and generates publication-ready visualizations.

---

## Method Overview

The workflow consists of four major steps.

### 1. Extract active HOR methylation

Using the input bgzip-compressed CpG methylation file together with active HOR annotations, the pipeline efficiently extracts methylation records overlapping each active HOR region and its flanking sequences using **tabix**.

### 2. Construct methylation profiles

Extracted CpG methylation records are aggregated into fixed genomic bins (user-defined, default: **1 kb**) to generate continuous methylation profiles across each active HOR region.

### 3. Identify candidate CDRs

For each active HOR, the pipeline estimates the representative methylation level and identifies genomic bins exhibiting substantial methylation depletion based on a user-defined threshold.

### 4. Merge adjacent bins

Neighboring hypomethylated bins are merged into continuous candidate CDRs according to the maximum allowed gap size and the minimum number of consecutive bins.

The pipeline reports:

- Candidate CDR coordinates
- Binned methylation profiles
- Smoothed methylation profiles
- QC figures


---

## Workflow

```text
bgzip methylation
      + active HOR annotation
                │
                ▼
     Low-depth filtering
                │
                ▼
    Methylation profile
                │
                ▼
   Candidate CDR detection
                │
                ▼
      Adjacent bin merging
                │
                ▼
     Candidate CDRs + QC
```

---

## Dependencies

### Required R packages

- data.table
- ggplot2
- scales
- optparse
- zoo

---

## Usage

```bash
Rscript call_CDR_from_active_HOR.R \
    --sample PAN010_Mat \
    --meth_file binned_methylation.tsv \
    --active_hor active_hor_regions.bed \
    --out_dir output_directory \
    --bin_gap 50000 \
    --min_bins 5 \
    --cutoff_frac 0.15 \
    --smooth_k 10 \
    --depth_cutoff 10
```

### Optional parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--bin_gap` | 50000 | Maximum gap allowed when merging adjacent candidate bins |
| `--min_bins` | 5 | Minimum number of consecutive bins required to define a CDR |
| `--cutoff_frac` | 0.15 | Fraction of the active HOR median methylation used to identify candidate CDRs |
| `--smooth_k` | 10 | Smoothing window size for visualization |
| `--depth_cutoff` | 10 | Minimum sequencing depth per genomic bin |

---

## Output

For each sample, the pipeline generates:

- Candidate CDR coordinates (BED)
- Binned methylation profiles
- Smoothed methylation profiles
- QC figures (PDF)

---
