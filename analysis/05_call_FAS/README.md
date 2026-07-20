# Module 05: Calculate FAS (FIRE Area Score)

This module calculates the **FIRE Area Score (FAS)**, a quantitative metric
developed in this study to measure the magnitude and genomic extent of Fiber-seq
inferred regulatory element (FIRE) enrichment within centromeric dip regions
(CDRs).

FAS compares FIRE coverage inside a CDR with the surrounding non-CDR portion of
the same active higher-order repeat (HOR) array. It combines relative FIRE
enrichment with the number of covered base pairs, providing a cumulative measure
of centromeric chromatin accessibility.

## Method overview

For each annotated active-HOR/CDR interval, the workflow:

1. Extracts Fiber-seq coverage records overlapping the active HOR with `tabix`.
2. Divides the active HOR into fixed-width genomic bins (5 kb by default).
3. Uses interval-overlap lengths to calculate the FIRE, linker, and nucleosome
   contributions in each bin.
4. Calculates the FIRE fraction and smooths it with a centered rolling window.
5. Estimates the non-CDR background after Winsorization.
6. Winsorizes FIRE fractions inside the CDR and calculates FAS.

The implemented score is:

```text
FAS = sum(((CDR-bin FIRE fraction - background FIRE fraction) /
           (background FIRE fraction + 1e-6)) * covered bp)
```

The workflow also reports a signed log transformation,
`sign(FAS) * log2(abs(FAS) + 1)`, and the ratio of mean CDR FIRE to background
FIRE. The implementation in this module should be used to reproduce the values
reported in the accompanying study.

## Workflow

```text
bgzip/tabix-indexed Fiber-seq coverage tracks
                         |
                         v
            Extract active-HOR records
                         |
                         v
          Length-weighted genomic binning
                         |
                         v
          Smooth binned FIRE fractions
                         |
                         v
       Estimate Winsorized non-CDR background
                         |
                         v
                  Calculate FAS
                         |
                         v
        Results, statistics, and QC plots
```

## Requirements

### Command-line software

- Bash 4.2 or later
- [HTSlib](https://www.htslib.org/) (`tabix`)
- R 4.0 or later

### R packages

- `data.table`
- `ggplot2`
- `zoo`
- `scales`

Install the required R packages with:

```r
install.packages(c("data.table", "ggplot2", "zoo", "scales"))
```

## Input files

### 1. Active-HOR/CDR annotation

Provide a five-column, tab-delimited file:

```text
#chr  activeHOR_start  activeHOR_end  CDR_start  CDR_end
chr9  45212770         47953550       45900000   46120000
```

Coordinates use the BED convention: zero-based starts and half-open ends. Each
CDR must fall entirely within its corresponding active-HOR interval.

### 2. Fiber-seq coverage manifest

Provide a two-column, tab-delimited sample manifest:

```text
#sampleid       data_path
PAN027Pat_iPSC  example_data/PAN027Pat_iPSC_chr9_activeHOR_all_element_coverages.bed.gz
PAN027Pat_NPC   example_data/PAN027Pat_NPC_chr9_activeHOR_all_element_coverages.bed.gz
```

Sample IDs must be unique. Relative data paths are interpreted from the
directory where `call_FAS.sh` is launched.

Each coverage file must be coordinate-sorted, bgzip-compressed, and accompanied
by a `.tbi` or `.csi` tabix index. Its first six columns must be:

```text
chr  start  end  FIRE_coverage  linker_coverage  nucleosome_coverage
```

Comment lines beginning with `#` and blank lines are ignored in both metadata
files.

## Usage

Run the example from the module directory:

```bash
bash call_FAS.sh \
  --fiberseq-list example_data/fs_files.tsv \
  --centromere example_data/activeHOR_CDR.tsv \
  --output-dir example_output
```

The script contains Slurm directives and can also be submitted with `sbatch`:

```bash
sbatch call_FAS.sh \
  --fiberseq-list example_data/fs_files.tsv \
  --centromere example_data/activeHOR_CDR.tsv \
  --output-dir example_output
```

Run `bash call_FAS.sh --help` for complete command-line help.

## Parameters

The descriptive long options and original version 1.0 options are both
supported. Defaults reproduce the version 1.0 pipeline.

| Long option | Original option | Default | Description |
|---|---|---:|---|
| `--output-dir` | `-o` | `call_FAS_output` | Output directory |
| `--bin-size` | `-bin` | `5000` | Genomic bin size in bp |
| `--smooth-window` | `-s` | `10` | Centered smoothing window in bins |
| `--cdr-expansion` | `-exp` | `0` | Symmetric CDR expansion in bp; negative values contract the CDR |
| `--outside-cdr-cutoff` | `-out_CDR_cutoff` | `0.10` | Winsorization fraction outside the CDR |
| `--inside-cdr-cutoff` | `-in_CDR_cutoff` | `0.05` | Winsorization fraction inside the CDR |

Expanded CDR boundaries are restricted to the annotated active-HOR interval.
The legacy `-t`/`--threads` option remains accepted for command compatibility,
but the version 1.0 workflow is sequential and does not use this value.

## Output

```text
call_FAS_output/
|-- FAS_results.tsv
|-- FAS_bin_summary.tsv
|-- FAS_activeHOR_coverage.tsv
|-- call_FAS.log
|-- FAS_QC_plots/
|   |-- all_chr_FAS_QC.pdf
|   |-- FAS_cell_difference_stats.tsv
|   `-- FAS_cell_difference_summary.pdf
`-- intermediate_files/
    |-- active_HOR_CDR.validated.tsv
    |-- fiberseq_manifest.validated.tsv
    `-- activeHOR_coverage/
```

- `FAS_results.tsv` contains sample-level FAS, signed `log2_FAS`, FIRE
  enrichment, background and CDR FIRE values, covered base pairs, bin counts,
  and analysis parameters.
- `FAS_bin_summary.tsv` contains the length-weighted and smoothed bin-level
  profiles used in the calculation.
- `FAS_activeHOR_coverage.tsv` contains the extracted coverage records supplied
  to the R analysis.
- `FAS_QC_plots/` contains per-sample profiles and optional iPSC-versus-NPC
  summary statistics.

## Reproducibility notes

- Cell types are inferred from sample names containing `PBMC`, `iPSC`, `NPC`,
  `Monocyte`, or `Macrophage`.
- The paired iPSC-versus-NPC statistics use the sample-name prefix before the
  first underscore as the individual/assembly identifier.
- Reusing an output directory replaces files with matching names. Use a separate
  directory for every parameter set.

