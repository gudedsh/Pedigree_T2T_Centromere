# Calling candidate centromeric dip regions

`call_CDR_candidate.sh` identifies candidate centromeric dip regions (CDRs)
from per-site methylation data within annotated active higher-order repeat (HOR)
regions. It bins methylation measurements with BEDTools and passes the resulting
per-bin table to `r_code/call_CDR_from_active_HOR.R` for CDR calling.

## Requirements

- Bash 4.2 or later
- BEDTools
- HTSlib (`tabix` and `bgzip`)
- R (4.1 or later) with `data.table`, `ggplot2`, and `optparse`

Install or verify the R dependencies with:

```bash
Rscript -e 'install.packages(c("data.table", "ggplot2", "optparse", "getopt"))'
Rscript -e 'invisible(lapply(c("data.table", "ggplot2", "optparse", "getopt"), requireNamespace, quietly = FALSE))'
```

The script can be run interactively or submitted to Slurm. Adjust the embedded
`#SBATCH` resource settings to match the local cluster configuration.

## Input files

The active-HOR file is a tab-delimited BED file. Comment lines beginning with
`#` are allowed:

```text
#chr  activeHOR_start  activeHOR_end
chr9  45249191         47361694
```

The methylation manifest has four tab-delimited columns. Column numbers are
one-based and refer to columns in the indexed methylation BED file:

```text
#sampleID                meth_column  depth_column  data_path
PAN027_Mat_PBMC_Prim_PB  4            6             ./example_data/PAN027_Mat_chr9_activeHOR_meth.bed.gz
```

Each methylation file must be coordinate-sorted, bgzip-compressed, and tabix
indexed. For example:

```bash
bgzip methylation.bed
tabix -p bed methylation.bed.gz
```

Sample names must be unique and may contain letters, numbers, periods,
underscores, and hyphens. Relative data paths in the manifest are interpreted
relative to the directory where the calling script is launched.

## Quick test with the example data

Run this command from the directory containing the script:

```bash
bash call_CDR_candidate.sh \
  --active-hor example_data/PAN027_Mat_chr9_activeHOR_annotation.tsv \
  --methylation-list example_data/meth_file.tsv \
  --output-dir example_output
```

On a Slurm cluster, use:

```bash
sbatch call_CDR_candidate.sh \
  --active-hor example_data/PAN027_Mat_chr9_activeHOR_annotation.tsv \
  --methylation-list example_data/meth_file.tsv \
  --output-dir example_output
```

Use `bash call_CDR_candidate.sh --help` to see all parameters. Add
`--keep-intermediates` to retain the binned methylation table passed to R; this
is useful for troubleshooting and method validation.

## Notes for reproducible use

- The default bin size is 5 kb.
- A bin must contain at least three mapped methylation records by default.
- The default minimum CDR length is five qualifying bins.
- Existing output files produced by the R script may be overwritten if the same
  output directory and sample names are reused. Use a new output directory for
  each parameter set.
- CDRs are called independently within each active-HOR interval. A candidate bin
  must have methylation below `median * (1 - cutoff_frac)`; the default is 15%
  below the median methylation of its active-HOR interval.
- Methylation values are expected to be percentages for the QC plot. The calling
  rule itself is scale-independent as long as all values use the same scale.

