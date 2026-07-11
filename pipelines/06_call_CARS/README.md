# Module 06: Call CARS (CENP-A Associated Region States)

This directory contains the computational pipeline for calculating the **CARS (CENP-A Associated Region States)** score and tracking CENP-A deposition enrichment. Utilizing long-read **CENP-A DiMeLo-seq** single-molecule mapping data, this workflow extracts localized modification signals (e.g., $m^6A$ or $m^6A/\text{AT}$ ratio), performs standardized window smoothing, implements robust dual-zone outlier trimming, and quantifies integrated footprint enrichment across **centromere-to-arm boundaries** and active Higher-Order Repeats (HORs).

---

## Workflow & Algorithmic Core

The main pipeline wrapper (`call_CARS_pipeline_V1.0.sh`) orchestrates parallel single-molecule feature extractions via Fibertools (`ft extract`), channels raw alignments into standard genomic window blocks via `bedtools intersect`, and routes the synchronized matrix into the core R algorithm (`call_cars_mars.r`) which executes across four modular phases:

1. **Signal Conditioning via Rolling Average**: Evaluates flanking macro-epigenetic trends over chromosome sequences using `zoo::rollmean` with a moving step parameter controlled by `--smooth_k` to generate a robust baseline.
2. **Dual-Zone Asymmetric Winsorization**: Symmetrically compresses high-frequency coverage peaks or non-specific hyper-methylation outliers using independent sample quantiles to stabilize regional reference lines:
   * **Flanking Background Bins (Outside CDR)**: Outliers are capped using the `--out_CDR_cutoff` fraction (Default: top/bottom 10%).
   * **Target Core Bins (Inside CDR)**: Outliers are capped using the `--in_CDR_cutoff` fraction (Default: top/bottom 5%).
3. **CARS Score Formulations**: Establishes the robust non-CDR baseline background mean ($\text{Background}$) from the trimmed flanking array, and integrates the relative enrichment delta scaled directly against the physical nucleotide width ($\text{Bin\_Width}$) of each target window:
   $$\text{Delta} = \frac{\text{Signal}_{\text{smooth}} - \text{Background}}{\text{Background} + 10^{-6}}$$
   $$\text{CARS} = \sum (\text{Delta}_{\text{cdr\_bin}} \times \text{Bin\_Width})$$
   The pipeline records absolute integrated $\text{CARS}$ scores, log-scaled scaling profiles ($\text{sign(CARS)} \times \log_2(\vert{}\text{CARS}\vert{} + 1)$), and un-winsorized raw metrics for comprehensive downstream reporting.

---

## Dependencies

Ensure the following infrastructure tools and R statistical packages are active within your cluster node environment:
* **Fibertools-rs** (`ft` binary for single-molecule long-read modification extraction)
* **Samtools / Bedtools** (Coordinate sorting, selective interval partitioning, and grouping)
* **R Environment** (>= 4.0) with required libraries:
  * `data.table` (High-efficiency stream matrix reading and calculation)
  * `zoo` (Infrastructure rolling computations)
  * `ggplot2`, `scales`, & `patchwork` (Aesthetic vector-grade composite layout rendering)

---

## Usage

### Command-Line Execution

Run the complete pipeline directly by invoking the master wrapper to automate parsing, overlap-splitting, and mathematical scoring:

```bash
bash call_CARS_pipeline_V1.0.sh \
  -dml /path/to/dimelo_bam_manifest.tsv \
  -cen /path/to/active_hor_cdr.bed \
  -o /path/to/output_directory \
  --bin_size 5000 \
  --smooth_k 10 \
  --cars_col m6A_per_AT \
  --cdr_expansion 0 \
  --out_CDR_cutoff 0.10 \
  --in_CDR_cutoff 0.05
