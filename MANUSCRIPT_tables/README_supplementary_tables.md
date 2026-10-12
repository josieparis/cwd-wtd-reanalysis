# Supplementary Tables — README

Supplementary Tables 1–5 for *Confounding of relatedness and exposure invalidates current polygenic genomic
predictions for chronic wasting disease in white-tailed deer*, a re-analysis of Seabury et al. (2020, *G3*
10:1433) and Seabury et al. (2022, *G3* 12:jkac109).

Each table is provided as a CSV file (UTF-8, comma-separated, fields quoted). All five tables with their legends
are also collected in `supplementary_tables.md`.

## Files

| Table | File | Content | Built by |
|---|---|---|---|
| S1 | `tableS1_aphis_funding.csv` | USDA APHIS-VS cooperative-agreement funds awarded to CWD predictive genetics, 2021–2025 | `build_supp_tables.R` |
| S2 | `tableS2_cohorts.csv` | The 2020 and 2022 datasets: samples, phenotypes, covariates, markers, analysis samples | `build_supp_tables.R` |
| S3 | `tableS3_naive_pc.csv` | Genomic-control λ and associated-SNP counts for structure-naive and PC-adjusted scans | `build_supp_tables.R` |
| S4 | `tableS4_prediction_repro.csv` | Reproduction of the in-sample GBLUP accuracy in Seabury et al. (2020), Table 1 | `build_supp_tables.R` |
| S5 | `tableS5_gxemm_perfarm.csv` | Per-farm liability-scale heritability from the GxEMM FREE model, farm as environment | `build_table5.R` |

## Conventions

- **Datasets.** "2020" is the Seabury et al. (2020) deposit (Dryad doi:10.5061/dryad.xd2547dcw; 807 deer);
  "2022" is the Seabury et al. (2022) deposit (Dryad doi:10.5061/dryad.wh70rxwnt; 1,170 deer with CWD status).
  97% of the 2020 deer are also in the 2022 dataset (Supplementary Table 2), so the two are not independent.
- **Regions.** The region coded "East" in the deposits is labelled "Northeast" throughout, including farm names.
- **Prevalence** is the sample case fraction (cases / deer with CWD status). It is not a population prevalence.
- **±** is a standard error for heritabilities (S5) and a standard deviation across replicates for prediction
  metrics (S4).
- **P-values** are written in scientific notation (e.g. 2.0×10⁻¹²). Negative values use a true minus sign (−).
- **—** means not applicable.

## Table details

### Supplementary Table 1 — APHIS-VS predictive-genetics funding

Transcribed from the APHIS Veterinary Services table "CA Funds Awarded to Predictive Genetics".

| Column | Definition |
|---|---|
| Year | Award year. 2025* was awarded non-competitively by APHIS-VS. |
| Predictive-genetics funding | US dollars awarded to predictive-genetics projects that year |
| Total funding | Total cooperative-agreement (CA) funding that year, US dollars |
| Percent of total | Predictive-genetics funding / total funding |
| Organizations awarded | Number of state agriculture agencies (State Ag) and universities (Univ) awarded |

The Total row is the sum of the yearly rows. The source table gives the Total funding as $22,528,638.20, but its
yearly amounts sum to $22,582,638.20, which is the value reported here.

### Supplementary Table 2 — The 2020 and 2022 datasets

One row per attribute, one column per dataset. Computed from the deposited genotypes and metadata, except the
SNP totals, which are as reported by the original studies.

| Row | Definition |
|---|---|
| Deer with CWD status | Deer with a binary CWD diagnosis (one deposited 2022 deer has no status) |
| Binary CWD (cases / controls) | CWD-positive / non-detect deer |
| Ordinal tissue score | 0 = non-detect, 1 = lymph node positive, 2 = lymph node and obex positive (2020 only) |
| Farms | Farm identifiers deposited (2022 only) |
| Analysis samples | Deer used per analysis: all deer (no covariates); complete sex and age (region models); farms with ≥ 10 deer (farm models); farms with ≥ 10 deer and both cases and controls (GxEMM and LDAK-PCGC) |
| Autosomal SNPs: relationship matrix | SNPs in the PLINK 1.9 `--make-rel` GRM (chromosomes 1–29, `--geno 0.05`) |
| Autosomal SNPs: association scans | SNPs in the single-GRM association scans (Table 1); 2020 additionally filtered to minor-allele frequency ≥ 0.05 |
| Genotype-identical pairs | Pairs within a dataset with a GRM off-diagonal value ≥ 0.9 |
| Overlap with the other dataset | 2020 deer whose best genotype match in 2022 has identity-by-state (IBS) ≥ 0.90 (Fig. S1) |

### Supplementary Table 3 — Structure-naive and PC-adjusted association scans

A calibration reference for Table 1, which reports the single-GRM mixed models.

| Column | Definition |
|---|---|
| Model | Structure-naive: per-SNP least-squares regression with no relatedness term. + 10 genomic PCs: the same with the top 10 principal components of the genotypes (PLINK PCA) as covariates |
| λ | Genomic-control inflation factor: median association χ² (1 df) divided by its expected median (0.455) |
| SNPs P < 5×10⁻⁸ | SNPs at genome-wide significance |
| SNPs P < 5×10⁻⁵ | SNPs at the authors' polygenic threshold |
| PRNP P | P-value for the *PRNP* codon-96 marker (G96S_G_to_A_at286) |

SNP sets: the 120,113 autosomal SNPs of Table 1 (2020) and the 119,962 autosomal SNPs passing the 5%
missingness filter (2022).

### Supplementary Table 4 — Reproduction of Seabury et al. (2020), Table 1

GBLUP (rrBLUP `mixed.solve`) on the 2020 GRM with random k-fold cross-validation, 50 replicates. Fixed effects
are added back for held-out deer, and predictions ≥ 0.5 are classified as cases (the authors' rule). Each
configuration has two rows: this study (mean ± SD over replicates) and the original (as reported, four decimals).

| Column | Definition |
|---|---|
| Folds (k) | Number of cross-validation folds |
| Covariates | Fixed effects: none, or sex, age and region |
| AUC | Area under the ROC curve |
| Accuracy | Proportion of deer classified correctly |
| Sensitivity / Specificity | Proportion of cases / controls classified correctly |
| MCC | Matthews correlation coefficient |
| RMSE | Root-mean-square error of the predicted value against the observed 0/1 status |

### Supplementary Table 5 — Per-farm GxEMM heritability

Liability-scale heritability for each farm from the GxEMM FREE model with farm as the environment. The sample is
the 1,120 deer on the 14 farms with ≥ 10 deer and both cases and controls, with sex and age as covariates and the
PLINK GRM. Fitted with LDAK5 (PCGC, `binary = TRUE`). The observed-scale FREE model did not converge.

| Column | Definition |
|---|---|
| Farm, Region | Farm identifier (region and farm number) and its region |
| n | Deer from the farm in the analysis |
| Prevalence | Case fraction on the farm |
| Liability h² ± SE (K = …) | Farm-specific heritability and its standard error, at assumed population prevalence K = 0.05, 0.10 or 0.26 |

Values outside 0 to 1 are inadmissible. Most estimates are inadmissible or have a standard error larger than the
estimate: the per-farm genetic variances are not identifiable (Table 5).

## Reproducing the tables

Run from the project root (or set `PROJECT_ROOT`):

```bash
Rscript reanalysis2022/scripts/build_supp_tables.R   # Supplementary Tables 1-4
Rscript reanalysis2022/scripts/build_table5.R        # Table 5 and Supplementary Table 5
```

The scripts read the Dryad deposits and these intermediate outputs:

| Table | Inputs | Produced by |
|---|---|---|
| S2 | `samples.tsv`; `cov_1170.csv`; `.bim` files; `grm.rel.gz`, `grm2022.rel.gz`; `mlma_bin.mlma`; `figS1_max_per_animal.csv` | `00_build_dataset.py`, `build_bed_2022.py`, PLINK 1.9 `--make-rel`, GCTA `--mlma`, `figS1_concordance_compute.py` |
| S3 | `gwas_binary.tsv.gz`, `gwas_score.tsv.gz` (2020); `gwas2022_binary.tsv.gz` (2022) | `04_gwas_lambda.py`, `make_fig2_data.py` |
| S4 | `table1_reproduction.csv` | `12_table1_reproduce.R` |
| S5 | `gxemm_farm_results.rds` | `run_gxemm_farm_cluster.R` (LDAK5, on a Linux cluster), with the GRM from `build_K_from_plink.R` |

S1 is entered directly in `build_supp_tables.R` from the APHIS-VS table.
