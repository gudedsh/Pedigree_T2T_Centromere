# Module 03: Call CDR (Centromere Dip Region) Candidates

This directory contains the pipeline for identifying **Centromere Dip Regions (CDRs)** candidates. The workflow integrates methylation frequencies with active Higher-Order Repeats (HORs), calculates localized methylation drops, and outputs CDR candidate regions alongside quality control visualizations.

---

## Workflow Overview

The core wrapper script coordinates data preparation and executes the underlying R tool `call_CDR_from_active_HOR.R` to process the data through four systematic stages:

1. **Filtering by Depth**: Filters out genomic positions that fall below the user-defined sequencing coverage threshold (`--depth_cutoff`).
2. **K-Smooth Trend Line Calculation**: Computes a moving average across chromosomes using `zoo::rollmean` to generate smooth trend lines exclusively for visualization.
3. **Adaptive Thresholding (Haplotype/Sample Aware)**: 
   * Dynamically calculates the `median` methylation level across the provided active HOR target region.
   * Defines a customized background drop cutoff based on a fraction of median
   * Flags any windows where the localized methylation rate drops below cutoff.
4. **Window Merging & Clustering**: Chains continuous low-methylation bins together into single candidate blocks, tolerating internal gaps up to `--bin_gap` and filtering for clusters containing at least `--min_bins`.

---

## Dependencies

Ensure the following R packages are installed in your active environment:
* `data.table`
* `ggplot2`
* `scales`
* `optparse`
* `zoo`

---

## Usage

### Command-Line Execution

Run the calling module directly from your terminal by executing the R tool with explicit path arguments:

```bash
Rscript call_CDR_from_active_HOR.R \
  --sample PAN010_Mat \
  --meth_file /path/to/binned_methylation.tsv \
  --active_hor /path/to/active_hor_regions.bed \
  --out_dir /path/to/output_directory \
  --bin_gap 50000 \
  --min_bins 5 \
  --cutoff_frac 0.15 \
  --smooth_k 10 \
  --depth_cutoff 10
