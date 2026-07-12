# Module 02: FIRE Call Pipeline

This module provides a wrapper for running the **FIRE (Fiber-seq Inferred Regulatory Elements)** pipeline to identify chromatin accessibility from PacBio Fiber-seq data.

The wrapper automates execution of the official FIRE workflow, manages sample-specific outputs, and generates the FIRE coverage files required for downstream analyses, including **FAS (FIRE Area Score)** calculation.

---

## Method Overview

The workflow consists of three major steps.

### 1. Run the FIRE pipeline

The input aligned Fiber-seq BAM file is processed using the official FIRE pipeline to identify accessible chromatin footprints and classify genomic elements.

### 2. Generate FIRE coverage tracks

The pipeline produces genome-wide FIRE coverage profiles (`all_element_coverages.bed.gz`), which summarize the fraction of FIRE, linker, and nucleosome footprints across the genome.

### 3. Prepare downstream analyses

Output files are automatically organized for downstream analyses, including FAS calculation and visualization of chromatin accessibility across centromeric regions.

---

## Dependencies

### External software

- FIRE Pipeline (official Snakemake workflow)
- Samtools (≥1.9)

For installation instructions, see the official FIRE documentation:

https://fiberseq.github.io/fire/run.html

---

## Usage

```bash
bash run_FIRE_pipeline_v1.0.sh \
    --sample sample_name \
    --bam input_alignment.bam \
    --ref reference.fa \
    --outdir output_directory \
    --threads 10 \
    --chrsize chromosome_sizes.txt
```

### Parameters

| Parameter | Description |
|-----------|-------------|
| `--sample` | Sample name |
| `--bam` | Input aligned Fiber-seq BAM file |
| `--ref` | Reference genome FASTA |
| `--outdir` | Output directory |
| `--threads` | Number of CPU threads |
| `--chrsize` | Chromosome size file (optional) |

---

## Output

The pipeline generates:

- FIRE footprint annotations
- Genome-wide FIRE coverage tracks (`all_element_coverages.bed.gz`)
- Other standard outputs produced by the official FIRE pipeline

These files serve as the input for **Module 05: Call FAS (FIRE Area Score)**.

---
