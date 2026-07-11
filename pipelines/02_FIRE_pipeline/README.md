# Module 02: FIRE (Fiber-seq Integrative Regulatory Elements) Call Pipeline

This directory contains the pipeline wrapper for running **FIRE (Fiber-seq Integrative Regulatory Elements)** chromatin accessibility calling on long-read PacBio platforms. It handles reference and alignment indexing, automatically formats dynamic Snakemake pipeline YAML configurations, tracks resource load management, and reformats output BigWig (trackHub) files to line up with localized chromosome size definitions.

---

## Workflow Overview

The core runner script `run_FIRE_pipeline_v1.0.sh` executes through five systematic phases:
1. **POSIX Command Argument Parsing**: Decodes critical genomic paths and allows flexible execution variables with robust built-in fallbacks.
2. **Indexing Verifications**: Confirms that the designated Reference and BAM input tracks are globally indexed via `samtools`, preventing structural Snakemake step crashes.
3. **Workspace Initialization**: Creates standalone operational workspaces avoiding complex file movements (`mv`) that corrupt background input directories.
4. **On-the-Fly Dynamic Config Generation**: Automates the rendering of tracking manifests (`.tbl`) and parameter configuration scripts (`.yaml`) tailored to the run specifications.
5. **Downstream TrackHub BigWig Reformatting**: Converts and reformats raw output bigwig alignments into standardized outputs matching your local chromosome size maps (`_chrsize.txt`).

---

## Dependencies

Ensure the following tools and environment configurations are exported within your active cluster session:
* **FIRE Pipeline executable** (Snakemake workflow framework)
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
