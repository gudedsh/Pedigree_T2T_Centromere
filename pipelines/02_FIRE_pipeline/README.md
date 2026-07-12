# Module 02: FIRE Call Pipeline

This directory contains the pipeline wrapper for running **FIRE (Fiber-seq Inferred Regulatory Elements)** chromatin accessibility calling on long-read PacBio platforms. 

---

---

## Dependencies

Ensure the following tools and environment configurations are exported within your active cluster session:
* **FIRE Pipeline ** (Snakemake workflow framework from https://fiberseq.github.io/fire/run.html)
* **Samtools** (v1.9 or later)

---

## Usage

### Command-Line Parameter Layout

```bash
bash run_FIRE_pipeline_v1.0.sh \
  --sample <sample_prefix> \
  --bam <input_alignment.bam> \
  --ref <reference_assembly.fa> \
  --outdir <output_directory> \
  --threads [10] \
  --chrsize [optional_chrom_size.txt]
