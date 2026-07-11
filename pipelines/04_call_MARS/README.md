# Module 04: Call MARS (Methylation-Based Active Region States)

This directory contains the computational pipeline for calculating **MARS (Methylation-Based Active Region States)** scores and mapping single-molecule epigenetic boundary states. Optimized for long-read **Fiber-seq** and **DiMeLo-seq** datasets, the workflow tracks DNA methylation drops across critical structural landmarks—such as **centromere-to-arm boundaries** and active Higher-Order Repeats (HORs)—by combining `tabix` localized data extraction, spatial bin mapping, and a specialized R mathematical framework (`call_mars.r`).

---

## Workflow & Algorithmic Core

The wrapper pipeline (`call_MARS_v1.0.sh`) coordinates initial workspace structures and routes parsed genomic matrices into the core R engine (`call_mars.r`), which executes through four automated stages:

1. **Cell-Type Aware Sorting**: Automatically group and structure incoming tracks based on parsed lineage substrings (prioritizing `PBMC` ➔ `iPSC` ➔ `NPC` ➔ `Monocyte` ➔ `Macrophage` order) to ensure coherent comparative display indexing.
2. **K-Smooth Trend Line Interpolation**: Computes a continuous localized moving average over chromosome blocks via `zoo::rollmean` to capture macro-level methylation trends used exclusively for single-point tracking visualization.
3. **Dual-Zone Robust Winsorization**: To prevent localized coverage drops or hyper-methylated outliers from destabilizing the baseline, the script applies an asymmetric quantile trimming cutoff to independent functional blocks:
   * **Background Regions (Outside CDR)**: Symmetrically trimmed using the `--out_CDR_cutoff` fraction (Default: top/bottom 10%).
   * **Target Regions (Inside CDR)**: Symmetrically trimmed using the `--in_CDR_cutoff` fraction (Default: top/bottom 5%).
4. **MARS Score Calculation**: Establishes the robust background baseline methylation level ($\text{Background}$) from the trimmed out-CDR windows, and integrates the single-molecule mathematical delta across the entire inner CDR span:
   $$\text{Delta} = \text{Background} - \text{Methylation}_{\text{winsorized}}$$
   $$\text{MARS} = \sum (\text{Delta}_{\text{cdr\_bins}})$$
   The pipeline reports both raw integrated absolute $\text{MARS}$ accumulations and normalized $\log_2(\text{MARS} + 10^{-6})$ scaling profiles.

---

## Dependencies

Ensure the following system binaries and language libraries are activated in your HPC execution environment:
* **Bedtools** (v2.26.0 or later)
* **Samtools / Tabix**
* **R Environment** (>= 4.0) with package extensions:
  * `data.table` (High-throughput file matrix streaming)
  * `ggplot2` & `scales` (Aesthetic publication vector rendering)
  * `zoo` (Rolling infrastructure computations)
  * `optparse` (POSIX command argument parsing)

---

## Usage

### Command-Line Direct Invocation

Execute the complete integrated wrapper to automate target expansions, bin-tiling, tabix sorting, and R scoring scripts in one step:

```bash
bash call_MARS_v1.0.sh \
  -ref /path/to/reference.fa \
  -cen /path/to/active_hor_cdr.bed \
  -me /path/to/methylation_manifest.tsv \
  -o /path/to/output_directory \
  -bin 1000 \
  -n 3 \
  -dl 5 \
  -dh 100 \
  -s 10 \
  -exp 0 \
  -out_CDR_cutoff 0.10 \
  -in_CDR_cutoff 0.05
