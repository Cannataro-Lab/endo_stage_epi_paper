# Summers et al. — Stepwise Evolution of Driver Mutations from Endometrial Hyperplasia to Carcinoma

**Authors:** Mary F. Summers, Kira A. Glasmacher, Sem Asmelash, J. Nic Fisk, Jeffrey D. Mandell, Jeffrey P. Townsend, Vincent L. Cannataro

------------------------------------------------------------------------

## Requirements

All analyses and figures were performed using R 4.4.2 with the following packages:

| Package           | Version  |
|-------------------|----------|
| cancereffectsizeR | 2.10.2   |
| ces.refset.hg19   | 1.1.3    |
| bbmle             | 1.0.25.1 |
| data.table        | 1.18.0   |
| tidyverse         | 2.0.0    |
| ggrepel           | 0.9.6    |
| patchwork         | 1.3.2    |
| scales            | 1.4.0    |

Install cancereffectsizeR and its reference set from GitHub:

``` r
remotes::install_github("Townsend-Lab-Yale/cancereffectsizeR@v2.10.2")
remotes::install_github("Townsend-Lab-Yale/ces.refset.hg19@v1.1.3")
```

------------------------------------------------------------------------

## Input data

| File | Source | Notes |
|---------------------|---------------------------|------------------------|
| `VCF1_r.csv` | Li et al. 2021 (J Pathol 253:119–128) | **Restricted** — contact Li et al. |
| `VCF2_r.csv` | Li et al. 2021 (J Pathol 253:119–128) | **Restricted** — contact Li et al. |
| `TCGA_ucec_data.maf` | TCGA via GDC | Downloaded automatically by `get_TCGA_project_MAF()` if absent |
| `cptac.maf` | CPTAC-3 via GDC | Downloaded automatically by `get_TCGA_project_MAF()` if absent |
| `clinical.tsv` | GDC data portal (TCGA-UCEC) | Clinical data for TCGA samples |
| `gdc_manifest.2023-02-28.txt` | GDC data portal (CPTAC-3) | Manifest selecting the 102 CPTAC endometrial samples |
| `hg38ToHg19.over.chain` | UCSC Genome Browser | Liftover chain file for hg38→hg19 coordinate conversion |

------------------------------------------------------------------------

## How to run

1.  Install all required packages (see Requirements above).
2.  Place all required input files in `input_data/` (see Input data above).
3.  Source `analysis.R`.

Figures will be added to the `figures/` folder.
