# Pedigree_T2T_Centromere

Code and analysis pipelines accompanying our study on centromere epigenetic remodeling and genetic stability using pedigree telomere-to-telomere (T2T) genome assemblies.

This repository provides reproducible analysis pipelines, genome annotations, interactive browser resources, and public datasets for studying DNA methylation, chromatin accessibility, CENP-A occupancy, and de novo mutations across human centromeres.

---

## Repository Structure

```text
Pedigree_T2T_Centromere/
├── annotations/                  Genome annotations
├── data/                         Public datasets and download links
├── interactive_browser_resource/ WashU Epigenome Browser resources
└── pipelines/
    ├── 01_process_PacBio_raw_data/
    ├── 02_FIRE_pipeline/
    ├── 03_call_CDR_candidate/
    ├── 04_call_MARS/
    ├── 05_call_FAS/
    ├── 06_call_CARS/
    └── 07_call_de_novo_mutation/
```

---

## Analysis Pipelines

| Module | Description |
|---------|-------------|
| **01** | Process raw PacBio long-read data, including Fiber-seq and DiMeLo-seq. |
| **02** | Run the FIRE pipeline to identify Fiber-seq inferred regulatory elements. |
| **03** | Identify Centromere Dip Region (CDR) candidates. |
| **04** | Calculate MARS (Methylation Area Score). |
| **05** | Calculate FAS (FIRE Area Score). |
| **06** | Calculate CARS (CENP-A Area Score). |
| **07** | Identify high-confidence de novo and somatic mutations. |

---

## Public Resources

### Genome assemblies

Phased T2T genome assemblies are publicly available through Zenodo.

**DOI**

https://doi.org/10.5281/zenodo.21326144

### Genome annotations

Assembly-specific annotations are available in the `annotations/` directory, including:

- Active HOR and CDR annotations
- CenSat annotations
- Coding exon annotations
- Promoter annotations

### WashU Epigenome Browser

Interactive visualization resources and tutorials are provided in:

```text
interactive_browser_resource/
```

---

## Citation

If you use this repository or associated datasets in your research, please cite:

> Dong S. *et al.* Fully T2T pedigree assemblies reveal genetic stability and epigenetic plasticity of human centromeres across inheritance and cell-fate transitions.

---

## License

The source code in this repository is distributed under the MIT License.

Genome assemblies and other datasets are distributed separately under the licenses specified in their respective repositories (e.g., Zenodo).

---


## Contact

**Shihua Dong**

dongsh0101@gmail.com

Department of Genetics

Washington University School of Medicine

St. Louis, Missouri, USA
