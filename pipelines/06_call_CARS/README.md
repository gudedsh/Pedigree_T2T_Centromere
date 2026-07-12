# Module 06: Call CARS (CENP-A Area Score)

This module implements **CARS (CENP-A Area Score)**, a quantitative metric developed in our study to measure CENP-A enrichment within centromere dip regions (CDRs).

Rather than relying solely on average CENP-A signal, CARS integrates both the **magnitude** and **genomic extent** of CENP-A enrichment, providing a robust quantitative measurement of centromeric CENP-A occupancy across different cell types and experimental conditions.

The pipeline accepts long-read CENP-A DiMeLo-seq BAM files together with CDR and active HOR annotations. It automatically extracts active HOR regions, quantifies CENP-A signals, calculates CARS scores, performs quality control, and generates publication-ready visualizations.

In addition to CARS calculation, the pipeline exports paired single-molecule CENP-A and CpG methylation measurements for downstream integrative analyses.

---

## Method Overview

The workflow consists of five major steps.

### 1. Extract active HOR regions

Reads overlapping active HOR regions and their flanking sequences are extracted from aligned DiMeLo-seq BAM files.

### 2. Quantify CENP-A signals

Single-molecule CENP-A signals are extracted using **Fibertools (`ft extract`)**. Signal metrics (e.g., `m6A_per_AT`) are summarized within fixed genomic bins across each active HOR region.

### 3. Extract single-molecule CpG methylation

In addition to CENP-A signals, the pipeline extracts CpG methylation levels for every genomic bin from each individual read. These single-molecule methylation profiles are generated simultaneously with CENP-A quantification, enabling downstream integrative analyses of CpG methylation and CENP-A occupancy at single-molecule resolution.

### 4. Estimate background signal

Background CENP-A signal is estimated from regions outside the CDR after robust Winsorization, reducing the influence of local outliers and sequencing coverage fluctuations.

### 5. Calculate CARS

For each CDR, the pipeline estimates a robust background CENP-A signal from the surrounding non-CDR regions and quantifies the cumulative CENP-A enrichment across the entire CDR relative to this background.

The pipeline reports:

- Raw CARS
- log₂-transformed CARS
- Background CENP-A signal
- CDR CENP-A signal
- CENP-A enrichment
- QC figures

---

## Workflow

```text
DiMeLo-seq BAM
      + active HOR annotation
                │
                ▼
      Extract active HOR reads
                │
                ▼
    ft extract (CENP-A signals)
                │
                ├──────────────► Single-read CpG methylation
                │                     │
                │                     ▼
                │         Single-molecule CENP-A /
                │         methylation analyses
                │
                ▼
      Fixed-bin CENP-A profile
                │
                ▼
 Background signal estimation
      (Robust Winsorization)
                │
                ▼
         CARS calculation
                │
                ▼
      QC plots & summary tables
```

---

## Dependencies

### External software

- Fibertools-rs
- Samtools
- Bedtools
- R (≥4.0)

### Required R packages

- data.table
- ggplot2
- zoo
- patchwork
- scales

---

## Usage

```bash
bash call_CARS_pipeline_V1.0.sh \
    -dml dimelo_manifest.tsv \
    -cen active_hor_cdr.bed \
    -o output_directory
```

### Optional parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--bin_size` | 5000 | Bin size (bp) |
| `--smooth_k` | 10 | Smoothing window size |
| `--cars_col` | `m6A_per_AT` | Signal used for CARS calculation |
| `--cdr_expansion` | 0 | CDR expansion size (bp) |
| `--out_CDR_cutoff` | 0.10 | Winsorization cutoff outside the CDR |
| `--in_CDR_cutoff` | 0.05 | Winsorization cutoff inside the CDR |

---

## Output

For each sample, the pipeline generates:

- Processed CENP-A signal profiles
- Single-molecule CpG methylation profiles
- CDR CENP-A enrichment profiles
- Background CENP-A estimates
- CARS summary tables
- QC figures (PDF)

---
