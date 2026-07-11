# Pedigree_T2T_Centromere

This repository contains the integrated bioinformatics pipeline for processing raw PacBio long-read sequencing data, specifically optimized for **Fiber-seq** and **DiMeLo-seq** datasets. The workflow features haplotype-aware alignment, competitive read phasing, and single-molecule methylation calling, making it highly suitable for resolving complex genomic regions such as **centromere-to-arm boundaries** in T2T (Telomere-to-Telomere) assemblies across pedigree lineages.

---

## Pipeline Workflow

The main wrapper script `run_process_PacBio_raw_data_pipeline_V1.0.sh` automates the execution of four sequential modules:

1. **Step 1: Predict m6A and QC** (`01_predict_m6a_qc.sh`)
   * Takes raw PacBio BAM files to predict $m^6A$ kinetics and outputs intermediate BAM files with modification tags alongside quality control metrics.
2. **Step 2: Align Unphased Reads to Haplotypes** (`02_align_unphased_to_haps.sh`)
   * Coordinates parallel alignment of the unphased kinetic reads to both maternal (`Mat`) and paternal (`Pat`) reference genomes.
3. **Step 3: Extract Scores and Phase Reads** (`03_phase_reads.sh`)
   * Extracts competitive alignment metrics and modification scores to phase individual long reads into distinct parental haplotype buckets based on user-defined thresholds.
4. **Step 4: Realign Phased Reads and Call Methylation** (`04_realign_phased_call_meth.sh`)
   * Performs a final, high-accuracy realignment of the partitioned maternal and paternal reads back to their respective native references and executes resolution-critical methylation calling.

---

## Directory Structure of Outputs

Upon successful completion, the pipeline automatically organizes outputs into the following structure within your designated `--outdir`:

```text
output_dir/
├── 01_predict_m6a_qc/
│   └── [sample].predict_m6a.bam
├── 02_align_unphased/
│   ├── Mat/
│   │   └── [sample].Mat.unphased.sorted.bam
│   └── Pat/
│       └── [sample].Pat.unphased.sorted.bam
├── 03_phase_reads/
│   ├── [sample].MatReads.bam
│   └── [sample].PatReads.bam
└── 04_realign_phased_call_meth/
    ├── Mat/
    │   └── [sample]_Mat.sorted.bam  <- Final Maternal Phased & Called Methylation
    └── Pat/
        └── [sample]_Pat.sorted.bam  <- Final Paternal Phased & Called Methylation

```


# Check help information
bash run_process_PacBio_raw_data_pipeline_V1.0.sh -h

# Full workflow execution example
bash run_process_PacBio_raw_data_pipeline_V1.0.sh \
  --input-bam /path/to/raw_input.bam \
  --sample PAN010 \
  --mat-ref /path/to/maternal_reference.fa \
  --pat-ref /path/to/paternal_reference.fa \
  --outdir /path/to/output_directory \
  --threads 16 \
  --min-len 1000 \
  --mg-cutoff 0


