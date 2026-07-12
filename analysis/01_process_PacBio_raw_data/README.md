# Module 01: Process PacBio Raw Data

This module provides the primary preprocessing pipeline for **PacBio HiFi**, **Fiber-seq**, and **DiMeLo-seq** datasets.

Starting from raw PacBio BAM files, the pipeline performs optional modification prediction, haplotype-aware alignment, competitive read phasing, and haplotype-specific methylation calling, producing phased BAM files and CpG methylation tracks for downstream centromere analyses.

---

## Method Overview

The workflow consists of four major steps.

### 1. Predict m6A and perform quality control *(optional)*

For Fiber-seq and DiMeLo-seq datasets, the pipeline predicts **m6A modifications** and nucleosome footprints from raw PacBio reads while generating quality control metrics.

This step can be skipped for standard PacBio HiFi whole-genome sequencing data that do not contain m6A modification signals.

### 2. Align unphased reads to both haplotypes

Raw reads are independently aligned to the maternal and paternal reference assemblies to generate competitive alignment candidates for each read.

### 3. Phase reads

Competitive alignment scores are evaluated to assign each read to its most likely parental haplotype, producing haplotype-partitioned read sets.

### 4. Realign phased reads and call methylation

Maternal and paternal read sets are realigned to their corresponding reference assemblies, followed by CpG methylation calling to generate final haplotype-resolved methylation tracks.

---

## Workflow

```text
Raw PacBio BAM
        │
        ▼
Predict m6A & QC (optional)
        │
        ▼
 Align to maternal & paternal
        │
        ▼
      Read phasing
        │
        ▼
 Realign phased reads
        │
        ▼
 CpG methylation calling
        │
        ▼
Phased BAM + methylation tracks
```

---

## Output Structure

```text
output_dir/
├── 01_predict_m6a_qc/
│   └── sample.predict_m6a.bam
│
├── 02_align_unphased/
│   ├── Mat/
│   └── Pat/
│
├── 03_phase_reads/
│   ├── sample.MatReads.bam
│   └── sample.PatReads.bam
│
└── 04_realign_phased_call_meth/
    ├── Mat/
    │   ├── sample_Mat.sorted.bam
    │   └── sample_Mat.methylation.bed.gz
    │
    └── Pat/
        ├── sample_Pat.sorted.bam
        └── sample_Pat.methylation.bed.gz
```

---

## Usage

### Show help

```bash
bash run_process_PacBio_raw_data_pipeline_V1.0.sh -h
```

### Run the complete pipeline

```bash
bash run_process_PacBio_raw_data_pipeline_V1.0.sh \
    --input-bam raw_input.bam \
    --sample PAN010 \
    --mat-ref maternal_reference.fa \
    --pat-ref paternal_reference.fa \
    --outdir output_directory \
    --threads 8 \
    --min-len 1000 \
    --mg-cutoff 0
```

### Main parameters

| Parameter | Description |
|-----------|-------------|
| `--input-bam` | Input PacBio BAM file |
| `--sample` | Sample name |
| `--mat-ref` | Maternal reference assembly |
| `--pat-ref` | Paternal reference assembly |
| `--outdir` | Output directory |
| `--threads` | Number of CPU threads |
| `--min-len` | Minimum read length |
| `--mg-cutoff` | Minimum alignment score difference used for read phasing |

---

## Output

For each sample, the pipeline generates:

- Haplotype-resolved BAM files
- Haplotype-resolved CpG methylation tracks
- m6A prediction results (Fiber-seq/DiMeLo-seq only)
- Quality control reports

The generated methylation tracks can be used directly as input for downstream analyses, including **CDR identification**, **MARS calculation**, and other centromere epigenomic analyses.
