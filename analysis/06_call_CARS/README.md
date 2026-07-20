# Module 06: Calculate CARS (CENP-A Area Score)

This module calculates the **CENP-A Area Score (CARS)**, a quantitative metric
developed in this study to measure the magnitude and genomic extent of CENP-A
enrichment within centromeric dip regions (CDRs).

The workflow processes long-read CENP-A DiMeLo-seq BAM files carrying Fiber-seq
tags. In addition to CARS, it exports bin-level CENP-A, CpG methylation,
methylspanning-patch (MSP), and nucleosome measurements for downstream
single-molecule analyses.

## Method overview

For every sample, the pipeline:

1. Extracts reads overlapping the annotated active higher-order repeat (HOR)
   regions.
2. Divides each active HOR into fixed-width bins (5 kb by default) and counts
   reference A, T, A+T, and CpG sites in each bin.
3. Retains read-bin combinations in which a read covers the complete bin.
4. Uses `ft extract` to obtain m6A, mCG, MSP, and nucleosome annotations.
5. Calculates bin-level `m6A_per_AT`, `mCG_per_CpG`, MSP fraction, and
   nucleosome fraction.
6. Smooths `m6A_per_AT` with a centered rolling window.
7. Retains genomic bins observed in all samples being compared for an active-HOR
   interval.
8. Estimates the non-CDR CENP-A background after Winsorization and calculates
   CARS across the CDR.

The implemented score is:

```text
CARS = sum(((CDR-bin signal - background signal) /
            (background signal + 1e-6)) * bin width)
```

The default signal is `m6A_per_AT`. The workflow also reports
`sign(CARS) * log2(abs(CARS) + 1)` and the CDR-to-background fold enrichment.
The implementation in this module should be used to reproduce the values
reported in the accompanying study.

## Workflow

```text
Indexed CENP-A DiMeLo-seq BAM files + reference FASTA
                            |
                            v
                 Extract active-HOR reads
                            |
                            v
             Build full-coverage read-bin pairs
                            |
                            v
           ft extract: m6A, mCG, MSP, nucleosomes
                            |
                 +----------+----------+
                 |                     |
                 v                     v
       Binned CENP-A profile    Single-molecule features
                 |
                 v
        Winsorized non-CDR background
                 |
                 v
              Calculate CARS
                 |
                 v
          Results tables and QC plots
```

## Requirements

### Command-line software

- Bash 4.2 or later
- [SAMtools](https://www.htslib.org/)
- [BEDTools](https://bedtools.readthedocs.io/)
- [fibertools-rs](https://github.com/fiberseq/fibertools-rs) (`ft`)
- GNU `sort`
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

The pipeline does not activate a hard-coded Conda environment. Activate an
environment containing these dependencies before running it.

## Input files

### 1. DiMeLo-seq BAM manifest

Provide a two-column, tab-delimited file:

```text
#sampleid       bam_path
PAN027Pat_iPSC  example_data/PAN027Pat_iPSC_chr9_activeHOR_CENPA_DiMeLo-seq_aligned.bam
PAN027Pat_NPC   example_data/PAN027Pat_NPC_chr9_activeHOR_CENPA_DiMeLo-seq_aligned.bam
```

Each BAM must be aligned to the supplied reference, contain the Fiber-seq tags
required by `ft extract`, and have a `.bai` or `.csi` index. Sample IDs must be
unique.

### 2. Active-HOR/CDR annotation

Provide a five-column, tab-delimited file:

```text
#chr  activeHOR_start  activeHOR_end  CDR_start  CDR_end
chr9  45212770         47953550       45900000   46120000
```

Coordinates use the BED convention: zero-based starts and half-open ends. Each
CDR must fall entirely within its corresponding active-HOR interval.

### 3. Reference genome

Provide the same reference FASTA used for BAM alignment. Uncompressed FASTA
files can be indexed in advance with:

```bash
samtools faidx reference.fa
```

The example reference is gzip-compressed. To avoid known `bedtools getfasta`
compatibility problems with compressed FASTA files, the pipeline automatically
decompresses it into `intermediate_files/reference_cache/`, creates a `.fai`
index for the cached copy, and uses that copy for reference-base counting. The
reference sequence and CARS calculation are unchanged.

Comment lines beginning with `#` and blank lines are ignored in the manifest and
annotation. Relative paths are interpreted from the directory where
`call_CARS.sh` is launched.

## Usage

Run from the module directory:

```bash
bash call_CARS.sh \
  --dimelo-list example_data/dml_files.tsv \
  --centromere example_data/activeHOR_CDR.tsv \
  --reference example_data/PAN027Pat_chr9.fa.gz \
  --output-dir example_output
```

The script contains Slurm directives and can also be submitted with `sbatch`:

```bash
sbatch call_CARS.sh \
  --dimelo-list example_data/dml_files.tsv \
  --centromere example_data/activeHOR_CDR.tsv \
  --reference example_data/PAN027Pat_chr9.fa.gz \
  --output-dir example_output
```

Run `bash call_CARS.sh --help` for complete command-line help.

## Parameters

Descriptive long options and the original version 1.0 options are both
supported. Defaults reproduce the version 1.0 workflow.

| Long option | Original option | Default | Description |
|---|---|---:|---|
| `--output-dir` | `-o` | `output_results` | Output directory |
| `--bin-size` | `-bin` | `5000` | Genomic bin size in bp |
| `--threads` | `-t` | `$SLURM_CPUS_PER_TASK` or `8` | Threads for SAMtools, fibertools, and sorting |
| `--smooth-window` | `-s` | `10` | Centered smoothing window in bins |
| `--cdr-expansion` | `-exp` | `0` | Symmetric CDR expansion in bp; negative values contract the CDR |
| `--outside-cdr-cutoff` | `-out_CDR_cutoff` | `0.10` | Winsorization fraction outside the CDR |
| `--inside-cdr-cutoff` | `-in_CDR_cutoff` | `0.05` | Winsorization fraction inside the CDR |

Expanded CDR boundaries are restricted to the active-HOR interval.

## Output

```text
output_results/
|-- call_CARS_results_from_CENPA_DiMeLo.tsv
|-- CARS_MARS_bin_summary.tsv
|-- call_CARS.log
|-- call_CARS_QC_plots/
|   |-- call_CARS_QC_plots.pdf
|   |-- summary_plot.pdf
|   |-- per_chr_plot.pdf
|   `-- results_summary.tsv
|-- cached_activeHOR_bam/
`-- intermediate_files/
```

- `call_CARS_results_from_CENPA_DiMeLo.tsv` contains CARS, signed `log2_CARS`,
  fold enrichment, raw and Winsorized background/CDR signals, bin counts, and
  analysis parameters.
- `CARS_MARS_bin_summary.tsv` contains the bin-level CENP-A and CpG methylation
  measurements used by the R analysis.
- Per-sample intermediate tables retain read-bin m6A, mCG, MSP, nucleosome, and
  reference-sequence counts.

## mCG quality-control rule

The version 1.0 mCG correction is retained unchanged. When the observed mCG
count exceeds the reference CpG count for a read-bin, the corrected count is
capped at the reference count. Excess of one or two calls is labeled
`EDGE_mCG_plus2`; larger excess is labeled `FAIL_mCG_excess`. Unaffected bins are
labeled `PASS`. Both raw and corrected counts are retained.

## Reproducibility notes

- CARS uses `m6A_per_AT` by default, exactly as in version 1.0.
- The bundled chromosome 9 example uses
  `example_data/PAN027Pat_chr9.fa.gz` and its
  `example_data/PAN027Pat_chr9.fa.gz.fai` index.
- Only reads spanning an entire bin contribute to that bin.
- Only bins observed in every compared sample for an active-HOR interval are
  retained before score calculation.
- Existing cached active-HOR BAMs are reused. Remove the relevant cached BAM and
  index before rerunning if the source BAM or active-HOR annotation has changed.
- Only summaries for samples in the current manifest are combined, preventing
  stale results from a previous run from entering the final table.

