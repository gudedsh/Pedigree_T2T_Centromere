# Module 05: Call FAS (FIRE Access Score)

This directory contains the computational pipeline for calculating the **FAS (FIRE Access Score)** and conducting downstream comparative statistical analyses across different cell lineage states. Utilizing single-molecule long-read chromatin accessibility profiles (e.g., from Fiber-seq datasets), this workflow extracts localized molecular fractions, performs length-weighted window binning, and implements robust asymmetric outlier trimming to model epigenetic shifts across **centromere-to-arm boundaries** and active Higher-Order Repeats (HORs).

---

## Algorithmic & Statistical Workflow

The main pipeline wrapper (`call_FAS_pipeline_V1.0.sh`) coordinates multi-threaded `tabix` selective data streaming, which is then routed directly into the core mathematical engine (`call_fas_from_fire_coverage.r`) to process data through five automated stages:

1. **Length-Weighted Tiling & Overlaps**: 
   * Dynamically partitions active HOR intervals into standardized window arrays based on `--bin_size`.
   * Employs `data.table::foverlaps` to evaluate intersecting data fragments and computes length-weighted average fractions for structural sub-components: `pfire` (FIRE fraction), `plinker` (linker fraction), and `pnuc` (nucleosome fraction).
2. **K-Smooth Vector Conditioning**: Evaluates moving trends over chromosome arrays via `zoo::rollmean` with a central alignment window controlled by `--smooth_k` to strip single-window background tracking noise.
3. **Dual-Zone Asymmetric Winsorization**: Symmetrically caps high- and low-end signal fluctuations using localized quantiles to stabilize the background profile baseline:
   * **Flanking Background Bins (Outside CDR)**: Outliers are compressed using the `--out_CDR_cutoff` fraction (Default: top/bottom 10%).
   * **Target Core Bins (Inside CDR)**: Outliers are compressed using the `--in_CDR_cutoff` fraction (Default: top/bottom 5%).
4. **FAS Metric Formulations**: Establishes the robust non-CDR background level ($\text{Background}$) from the trimmed flanking array and integrates the fractional delta scaled against the width of each active bin:
   $$\text{Delta} = \frac{\text{FIRE\_w} - \text{Background}}{\text{Background} + 10^{-6}}$$
   $$\text{FAS} = \sum (\text{Delta}_{\text{bin}} \times \text{Covered\_bp})$$
   The R script outputs both raw integrated absolute $\text{FAS}$ indexes, normalized $\text{sign(FAS)} \times \log_2(\vert{}\text{FAS}\vert{} + 1)$ scaling, and standard `FIRE_enrichment` ratios.
5. **Paired Lineage Statistics**: Automatically groups samples into functional lineages (prioritizing `PBMC`, `iPSC`, `NPC`, `Monocyte`, `Macrophage` strings). It performs a **Paired Wilcoxon Signed-Rank Test** explicitly contrasting target states (such as `NPC` vs `iPSC`) across shared assembly coordinates, outputting Benjamin-Hochberg (`BH`) adjusted **FDR** p-values.

---

## Dependencies

Ensure the following system tools and R scripting architectures are active in your cluster runtime path:
* **Samtools / Tabix** (High-throughput coordinate sorting and indexing)
* **R Environment** (>= 4.0) with required mathematical libraries:
  * `data.table` (High-efficiency genomic data manipulation)
  * `zoo` (Infrastructure rolling computations)
  * `ggplot2` & `scales` (Vector graphic rendering layouts)

---

## Usage

### Command-Line Execution

Run the complete pipeline directly by invoking the master execution wrapper script:

```bash
bash call_FAS_pipeline_V1.0.sh \
  -fa /path/to/coverage_manifest.tsv \
  -cen /path/to/active_hor_cdr.bed \
  -o /path/to/output_directory \
  --bin_size 2000 \
  --smooth_k 10 \
  --cdr_expansion 0 \
  --out_CDR_cutoff 0.10 \
  --in_CDR_cutoff 0.05
