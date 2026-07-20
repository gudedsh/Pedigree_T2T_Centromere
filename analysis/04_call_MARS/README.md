# Module 04: Calculate MARS (Methylation Area Score)

This module calculates the **Methylation Area Score (MARS)**, a quantitative
metric developed in this study to measure the magnitude and genomic extent of
DNA hypomethylation within centromeric dip regions (CDRs).

MARS summarizes the cumulative methylation depletion across a CDR relative to
the surrounding non-CDR portion of the same active higher-order repeat (HOR)
array. The workflow supports bgzip-compressed, tabix-indexed CpG methylation
tracks generated from PacBio HiFi, Fiber-seq, Oxford Nanopore, or other
platforms with compatible BED-like output.

## Method overview

For every annotated active-HOR/CDR interval, the workflow:

1. Extracts CpG records overlapping the active HOR using `tabix`.
2. Retains CpGs within the specified depth range.
3. Aggregates methylation into fixed-width genomic bins (1 kb by default).
4. Smooths each sample's binned methylation profile with a centered rolling
   window.
5. Retains genomic bins observed in all samples being compared for that
   active-HOR interval.
6. Estimates background methylation from bins outside the CDR after
   Winsorization.
7. Winsorizes methylation values inside the CDR and calculates MARS as the
   cumulative difference between background and CDR-bin methylation.

In simplified form:

```text
MARS = sum(background methylation - CDR-bin methylation)
```

Both raw MARS and `log2(MARS + 1e-6)` are reported. The implementation in this
module should be used to reproduce the values reported in the accompanying
study.

## Workflow

```text
bgzip/tabix-indexed CpG methylation tracks
                       |
                       v
          Extract CpGs within active HORs
                       |
                       v
             Depth filter and bin CpGs
                       |
                       v
          Smooth binned methylation profiles
                       |
                       v
        Retain bins shared across all samples
                       |
                       v
     Estimate Winsorized non-CDR background
                       |
                       v
               Calculate MARS
                       |
                       v
             Results table and QC plots
```

## Requirements

### Command-line software

- Bash 4.2 or later
- [BEDTools](https://bedtools.readthedocs.io/) 2.26 or later
- [HTSlib](https://www.htslib.org/) (`bgzip` and `tabix`)
- GNU `shuf`
- R 4.0 or later

### R packages

- `data.table`
- `ggplot2`
- `zoo`
- `scales`

The required R packages can be installed with:

```r
install.packages(c("data.table", "ggplot2", "zoo", "scales"))
```

## Input files

### 1. Active-HOR/CDR annotation

The annotation is a tab-delimited file containing five columns:

```text
#chr  activeHOR_start  activeHOR_end  CDR_start  CDR_end
chr9  45364117         47474748       45560000   45950000
```

Coordinates use the BED convention: zero-based starts and half-open ends. Each
CDR must fall entirely within its corresponding active-HOR interval.

### 2. Methylation manifest

The sample manifest is a four-column, tab-delimited file:

```text
#sampleid        meth_column  depth_column  data_path
PAN027Pat_iPSC   4            6             example_data/PAN027Pat_iPSC_chr9_activeHOR_meth.bed.gz
PAN027Pat_NPC    4            6             example_data/PAN027Pat_NPC_chr9_activeHOR_meth.bed.gz
PAN027Pat_PBMC   4            6             example_data/PAN027Pat_PBMC_chr9_activeHOR_meth.bed.gz
```

`meth_column` and `depth_column` are one-based column numbers in the associated
methylation file. Sample IDs must be unique.

Each methylation file must:

- contain chromosome, start, and end in its first three columns;
- be sorted by genomic coordinate;
- be compressed with `bgzip`; and
- have a `.tbi` or `.csi` tabix index.

For example:

```bash
bgzip methylation.bed
tabix -p bed methylation.bed.gz
```

Comment lines beginning with `#` and blank lines are ignored in both input
metadata files. Relative methylation paths are interpreted from the directory
where `call_MARS.sh` is launched.

## Usage

Run the example from the module directory:

```bash
bash call_MARS.sh \
  --centromere example_data/activeHOR_CDR.tsv \
  --methylation-list example_data/meth_files.tsv \
  --output-dir example_output
```

The script contains Slurm directives and can also be submitted with:

```bash
sbatch call_MARS.sh \
  --centromere example_data/activeHOR_CDR.tsv \
  --methylation-list example_data/meth_files.tsv \
  --output-dir example_output
```

Run `bash call_MARS.sh --help` for the complete command-line help.

## Parameters

Both the descriptive long options and the original short options are supported.
The defaults reproduce the version 1.0 workflow.

| Long option | Original option | Default | Description |
|---|---:|---:|---|
| `--output-dir` | `-o` | `output_results` | Output directory |
| `--bin-size` | `-bin` | `1000` | Genomic bin size in bp |
| `--min-cpg` | `-n` | `3` | Minimum CpGs required per bin |
| `--min-depth` | `-dl` | `5` | Minimum per-CpG depth |
| `--max-depth` | `-dh` | `100` | Maximum per-CpG depth |
| `--smooth-window` | `-s` | `10` | Centered smoothing window in bins |
| `--cdr-expansion` | `-exp` | `0` | Symmetric CDR boundary expansion in bp; negative values contract the CDR |
| `--outside-cdr-cutoff` | `-out_CDR_cutoff` | `0.10` | Winsorization fraction outside the CDR |
| `--inside-cdr-cutoff` | `-in_CDR_cutoff` | `0.05` | Winsorization fraction inside the CDR |

Expanded CDR boundaries are restricted to the annotated active-HOR interval.

## Output

The output directory contains:

```text
output_results/
|-- MARS_results.tsv
|-- call_MARS.log
|-- MARS_QC_plots/
|   |-- <chromosome>_MARS_QC1_all_samples.pdf
|   `-- <chromosome>_MARS_QC2_sample_facet.pdf
`-- intermediate_files/
    |-- combined_active_HOR.<bin_size>bp_methylation.tsv.gz
    |-- active_HOR.<bin_size>bp_bins.bed
    `-- per-sample intermediate files
```

`MARS_results.tsv` reports the annotation, MARS measurements, methylation
summaries, bin counts, and analysis parameters for every sample and CDR. The QC
plots show the smoothed methylation profiles, annotated CDR, background level,
and sample-level MARS values.

## Reproducibility notes

- As in version 1.0, up to 1,000 methylation values are sampled with `shuf` to
  determine whether each input uses a 0-1 or 0-100 methylation scale. The sampled
  mean is recorded in `call_MARS.log`.
- Bins not observed in every sample for an active-HOR interval are excluded
  before MARS calculation.
- Reusing an output directory replaces files with matching names. Use a separate
  directory for each parameter set.

