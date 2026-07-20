# Pedigree T2T Centromere

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21459252.svg)](https://doi.org/10.5281/zenodo.21459252)

Code, annotations, example datasets, and interactive resources accompanying our
study of centromere epigenetic remodeling and genetic stability using
haplotype-resolved, telomere-to-telomere (T2T) pedigree genome assemblies.

This repository provides reproducible workflows for processing PacBio
Fiber-seq and CENP-A DiMeLo-seq data, identifying centromeric dip regions
(CDRs), quantifying centromeric DNA methylation and chromatin features, and
detecting de novo and somatic mutations. Assembly-specific annotations and
small example datasets are included to facilitate testing.

## Repository overview

```text
Pedigree_T2T_Centromere/
├── analysis/
│   ├── 01_process_PacBio_raw_data/
│   ├── 02_FIRE_pipeline/
│   ├── 03_call_CDR_candidate/
│   ├── 04_call_MARS/
│   ├── 05_call_FAS/
│   ├── 06_call_CARS/
│   ├── 07_call_de_novo_mutation/
│   └── 08_extract_footprint_and_cpg_methylation/
├── annotations/
│   ├── activeHOR_CDR/
│   ├── cenSat/
│   ├── coding_exon/
│   └── promoter/
├── data/
└── washu_browser/
```

Each analysis module has its own README with software requirements, input
formats, parameters, example commands, and output descriptions. Modules 01,
03–06, and 08 include compact example data; tested example outputs are also
provided where applicable.

## Analysis modules

| Module | Workflow | Description |
|---:|---|---|
| 01 | [Process PacBio raw data](analysis/01_process_PacBio_raw_data/) | Process PacBio Fiber-seq and DiMeLo-seq reads, including m6A prediction, alignment, phasing, realignment, and methylation calling. |
| 02 | [FIRE pipeline](analysis/02_FIRE_pipeline/) | Identify Fiber-seq inferred regulatory elements (FIREs) and generate FIRE coverage profiles. |
| 03 | [Call CDR candidates](analysis/03_call_CDR_candidate/) | Identify candidate centromeric dip regions from active-HOR DNA methylation profiles. |
| 04 | [Calculate MARS](analysis/04_call_MARS/) | Calculate the Methylation Area Score, which summarizes the magnitude and genomic extent of CDR hypomethylation. |
| 05 | [Calculate FAS](analysis/05_call_FAS/) | Calculate the FIRE Area Score, which summarizes FIRE enrichment across CDRs. |
| 06 | [Calculate CARS](analysis/06_call_CARS/) | Calculate the CENP-A Area Score from CENP-A DiMeLo-seq data and export associated bin-level methylation measurements. |
| 07 | [Call de novo mutations](analysis/07_call_de_novo_mutation/) | Identify high-confidence germline de novo and somatic mutations using DeepVariant and GATK Mutect2 workflows. |
| 08 | [Extract footprints and CpG methylation](analysis/08_extract_footprint_and_cpg_methylation/) | Extract inferred nucleosome footprints, footprint sizes, reference CpG content, and read-matched CpG methylation calls from tagged BAM files. |

## Quantitative centromere scores

Three complementary scores are implemented in this repository:

| Score | Molecular feature | Interpretation |
|---|---|---|
| **MARS** | DNA methylation | Cumulative methylation depletion within a CDR relative to the surrounding non-CDR active HOR. |
| **FAS** | Fiber-seq FIRE coverage | Cumulative FIRE enrichment within a CDR relative to the surrounding non-CDR active HOR. |
| **CARS** | CENP-A DiMeLo-seq signal | Cumulative CENP-A enrichment within a CDR relative to the surrounding non-CDR active HOR. |

For reproducibility, use the implementations and defaults documented in the
corresponding module READMEs.

## Getting started

Clone or download the repository, enter an analysis module, and begin with its
example dataset. For example:

```bash
cd analysis/04_call_MARS

bash call_MARS_v1.0.sh \
  --centromere example_data/activeHOR_CDR.tsv \
  --methylation-list example_data/meth_files.tsv \
  --output-dir test_output
```

The workflows were designed for Linux and HPC environments. Several scripts
contain Slurm resource directives but can also be executed directly with Bash.
Software requirements differ by module and include combinations of SAMtools,
BEDTools, HTSlib, fibertools-rs, R, DeepVariant, and GATK. Consult the module
README before running a workflow.

## Genome annotations

Assembly-specific annotations are available under [annotations/](annotations/):

- [Active-HOR and CDR annotations](annotations/activeHOR_CDR/)
- [Centromeric satellite annotations](annotations/cenSat/)
- [Unique coding-exon annotations](annotations/coding_exon/)
- [Promoter annotations](annotations/promoter/)

Annotations are provided for the maternal and paternal haplotypes of PAN010,
PAN011, PAN027, and PAN028. These files use assembly-specific coordinates; use
each annotation only with its matching genome assembly and aligned data.

## Public data and interactive resources

### T2T genome assemblies

The phased T2T pedigree genome assemblies are available through Zenodo:

**[Zenodo DOI: 10.5281/zenodo.21326144](https://doi.org/10.5281/zenodo.21326144)**

Additional dataset descriptions and download information are provided in
[data/](data/).

### WashU Epigenome Browser

Interactive visualization resources and a PDF tutorial are available in
[washu_browser/](washu_browser/).

## Reproducibility notes

- Use the example datasets to verify software installation before processing
  complete datasets.
- Keep genome assemblies, annotations, BAM files, and tabix-indexed tracks in
  the same coordinate system.
- Record software and package versions when reproducing the analyses.
- Use a separate output directory for every sample set or parameter set to
  avoid mixing cached and newly generated intermediate files.

## Citation

If you use this repository, its annotations, or the associated datasets, please
cite:

> Dong S. *et al.* Fully T2T pedigree assemblies reveal genetic stability and
> epigenetic plasticity of human centromeres across inheritance and cell-fate
> transitions.

Please update this citation with the final journal reference when available.

## License

Source code in this repository is distributed under the MIT License. Genome
assemblies and other externally hosted datasets are distributed under the terms
specified by their respective repositories.

## Contact

**Shihua Dong**  
Department of Genetics  
Washington University School of Medicine  
St. Louis, Missouri, USA  
[dongsh0101@gmail.com](mailto:dongsh0101@gmail.com)

