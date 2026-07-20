# Module 08: Extract nucleosome footprints and CpG methylation

This module extracts Fiber-seq-inferred nucleosome footprints and CpG
methylation calls from an aligned BAM file within annotated active higher-order
repeat (HOR) regions.

For every inferred nucleosome on every read, the output reports its genomic
coordinates, footprint size, number of reference CpG dinucleotides, and number
of methylated CpG calls observed on the same read.

## Method overview

The workflow:

1. Extracts reads overlapping the supplied active-HOR intervals with SAMtools.
2. Uses `ft extract` to obtain CpG calls and inferred nucleosome footprints.
3. Converts the BED12 blocks produced by fibertools into individual BED6
   intervals.
4. Uses `bedtools nuc -pattern CG` to count reference CpGs within each
   nucleosome footprint.
5. Intersects footprints with CpG calls and counts only calls whose read ID
   matches the footprint read ID.
6. Reports footprints without a methylated CpG call with `mcg = 0`.

The extraction and counting definitions are unchanged from the original
analysis script.

## Requirements

- Bash 4.2 or later
- [SAMtools](https://www.htslib.org/)
- [BEDTools](https://bedtools.readthedocs.io/)
- [fibertools-rs](https://github.com/fiberseq/fibertools-rs) (`ft`)
- `bgzip`, `gzip`, GNU `sort`, and AWK

## Input files

### Tagged BAM

The BAM must:

- be aligned to the supplied reference;
- contain the Fiber-seq tags required by `ft extract`; and
- have a `.bai` or `.csi` index.

The bundled example uses:

```text
example_data/PAN027Pat_chr9_activeHOR_NPC_aligned.bam
example_data/PAN027Pat_chr9_activeHOR_NPC_aligned.bam.bai
```

### Active-HOR/CDR annotation

Provide a five-column, tab-delimited file:

```text
#chr  activeHOR_start  activeHOR_end  CDR_start  CDR_end
chr9  45212770         47953550       45900000   46120000
```

This module extracts all active-HOR intervals in the supplied file. CDR
coordinates are validated and retained as an intermediate BED file, but they do
not alter footprint extraction.

### Reference FASTA

Use the same reference sequence employed for BAM alignment. The example uses:

```text
example_data/PAN027Pat_chr9.fa.gz
```

Some BEDTools builds cannot read compressed FASTA files reliably. When the
reference ends in `.gz`, the script automatically creates and indexes an
uncompressed copy under `intermediate_files/reference_cache/`. This does not
change the reference sequence or counting method.

In the fibertools CLI used by this workflow, `ft extract --reference` is a
coordinate-output flag and does not accept a FASTA argument. The FASTA path is
used separately by `bedtools nuc` for reference-CpG counting.

## Usage

Run the example from the module directory:

```bash
bash extract_footprint_cpg_methylation.sh \
  --bam example_data/PAN027Pat_chr9_activeHOR_NPC_aligned.bam \
  --regions example_data/activeHOR_CDR.tsv \
  --reference example_data/PAN027Pat_chr9.fa.gz \
  --sample PAN027Pat_NPC \
  --output-dir example_output
```

The script contains Slurm resource directives and can also be submitted with
`sbatch` using the same arguments.

## Parameters

| Parameter | Default | Description |
|---|---:|---|
| `--bam` | required | Tagged, aligned BAM file |
| `--regions` | required | Active-HOR/CDR annotation |
| `--reference` | required | Matching reference FASTA |
| `--sample` | required | Sample identifier |
| `--assembly` | inferred | Assembly label; defaults to the sample name before its final underscore |
| `--cell` | inferred | Cell label; defaults to the sample name after its final underscore |
| `--output-dir` | `footprint_cpg_output` | Output directory |
| `--threads` | Slurm CPUs or `8` | SAMtools, fibertools, and sorting threads |
| `--mapq` | `1` | Minimum read mapping quality |
| `--keep-intermediates` | off | Retain extracted BAM and feature BED files |

For `PAN027Pat_NPC`, the inferred assembly is `PAN027Pat` and the inferred cell
label is `NPC`. Use explicit `--assembly` and `--cell` options if a sample name
does not follow this convention.

## Output

The primary output is:

```text
<output_dir>/<sample>.nucleosome_CpG_methylation.tsv.gz
```

It contains the following columns:

| Column | Description |
|---|---|
| `chr` | Reference sequence name |
| `start` | Zero-based footprint start |
| `end` | Half-open footprint end |
| `nucleosome_size_bp` | Footprint length (`end - start`) |
| `ncg` | Number of reference `CG` dinucleotides within the footprint |
| `mcg` | Number of extracted CpG methylation calls on the same read within the footprint |
| `asm` | Assembly label |
| `cell` | Cell label |
| `read_id` | Read identifier |

The output is bgzip-compressed. A run log and validated/intermediate files are
written under the selected output directory. Add `--keep-intermediates` when the
BED12/BED6 feature files or active-HOR BAM are needed for inspection.

## Reproducibility notes

- CpG calls are matched to nucleosome footprints by both genomic overlap and
  read ID.
- The default mapping-quality threshold remains `MAPQ >= 1`.
- A BED12 block with a non-positive reported length is represented as one base,
  matching the original script.
- No methylation fraction or additional coverage filter is applied.

