# Genome Annotations

This directory contains genome annotation files used throughout the analysis pipelines in this repository.

All annotations are generated on the pedigree T2T assemblies and provide standardized genomic features required for centromere epigenomic analyses.

---

## Available Annotations

| Directory | Description |
|----------|-------------|
| `activeHOR_CDR/` | Active HOR annotations together with corresponding Centromere Dip Region (CDR) annotations. |
| `cenSat/` | CenSat repeat annotations for centromeric satellite sequences. |
| `coding_exon/` | Protein-coding exon annotations. |
| `promoter/` | Promoter annotations derived from protein-coding gene models. |

---

## Usage

These annotation files are used by multiple analysis modules, including:

- **Module 03:** Call CDR candidates
- **Module 04:** Call MARS
- **Module 05:** Call FAS
- **Module 06:** Call CARS
- **Module 07:** De novo / Somatic mutation annotation

Each annotation file is assembly-specific and should be matched with the corresponding reference assembly used for analysis.

---

## File Naming

Annotation files follow the naming convention:

```text
<assembly>_<annotation>.bed.gz
```

For example:

```text
PAN010Mat_activeHOR_CDR_annotation.tsv.gz
PAN010Pat_cenSat.bed.gz
PAN027Mat_gene_annotation.bed.gz
```

where `<assembly>` corresponds to an individual haplotype assembly (e.g., `PAN010Mat`, `PAN010Pat`, `PAN027Mat`).

---

## Notes

All genomic coordinates are based on the corresponding pedigree T2T assemblies described in this study.
