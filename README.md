#### GitHub repo containing scripts pertaining to the article:

## Confounding of relatedness and exposure invalidates current polygenic genomic predictions for chronic wasting disease in white-tailed deer

#### Josephine R Paris, Aaron Westmoreland, Marty Kardos, W. Chris Funk, Jennifer L Malmberg

Archived under the Zenodo DOI: xxx

----------------------------------------------------------

We re-analysed the data from two studies of chronic wasting disease (CWD) in farmed white-tailed deer:

| Dataset | Deer | Study | Data |
| ------------- | ------------- | ------------- | ------------- |
| 2020 | 807 | Seabury et al. 2020, G3 10:1433 | https://doi.org/10.5061/dryad.xd2547dcw |
| 2022 | 1,170 | Seabury et al. 2022, G3 12:jkac109 | https://doi.org/10.5061/dryad.wh70rxwnt |

781 of the 807 deer from 2020 are also in the 2022 dataset.

Download both Dryad deposits into the main project folder (`Seabury_data`). The R scripts set the working directory at the top (`setwd("<dir>/Seabury_data")`), the shell scripts use `MASTER=<dir>/Seabury_data`, and the python scripts are run from the main project folder.

### Folders

| Folder | Contents |
| ------------- | ------------- |
| `reanalysis/scripts` | 2020 dataset: data set-up, plink, PCA, association, genomic prediction |
| `reanalysis/02_assoc` | 2020 GEMMA and GCTA association runs |
| `reanalysis2022/scripts` | 2022 dataset analyses, and the scripts for the figures and tables |
| `reanalysis2022/cluster_h2` | GCTA and GEMMA heritability (run on the cluster) |
| `reanalysis2022/cluster_gxemm`, `cluster_gxemm_2020` | GxEMM gene-by-environment models (run on the cluster) |
| `reanalysis2022/mlma_results` | 2022 GCTA association runs (run on the cluster) |
| `MANUSCRIPT_tables` | table outputs and the supplementary table README |

GCTA and LDAK (used by GxEMM) were run on a Linux cluster. GxEMM is from https://github.com/andywdahl/gxemm.

### Scripts for each table and figure

| Manuscript | Script |
| ------------- | ------------- |
| Fig 1 and Results section 1 | `coresidency.R`, `fig1_structure.R` |
| Table 1 | `build_table1.R` (association runs in `reanalysis/02_assoc`, `cluster_h2/gemma`, `mlma_results`) |
| Table 2, Fig 2a,c | `prediction_lofo.R`, `build_table2.R`, `fig2_prediction_collapse.R` |
| Fig 2b | `relatedness_binned_auc.R` |
| Fig 2d | `within_vs_across_region.R` |
| Fig 3 | `prnp_vs_gebv_threeway.R`, `fig3_prnp_vs_gebv.R` |
| AUC ceiling | `auc_ceiling_wray.R` |
| Table 3 | `build_table3.R`, `h2_with_se.R` |
| Table 4 | `build_table4.R` |
| Table 5, Supplementary Table 5 | `build_table5.R` (GxEMM runs in `cluster_gxemm`) |
| Fig 4 | `fig4_heritability.R` |
| Supplementary Tables 1-4 | `build_supp_tables.R` |
| Supplementary Fig 1 | `figS1_concordance_compute.py`, `fig_s1_concordance.R` |
| Supplementary Fig 2, Supplementary Table 3 | `fig2_calibration.R`, `make_fig2_data.py`, `reanalysis/scripts/04_gwas_lambda.py` |
| Supplementary Fig 3 | `arsb_perfarm.R`, `fig_s3_arsb.R` |
| Table S1 recount, threshold counts | `reanalysis/scripts/13_threshold_counts.py` |
| Tissue score as a covariate | `tissue_score_covariate.R` |
| Within-farm attenuation | `within_farm_attenuation.R` |
| Winner's curse, PRNP effect size | `reanalysis/scripts/06_winners_curse.py` |
| 2020 prediction (Supplementary Table 4) | `reanalysis/scripts/12_table1_reproduce.R` |

### Software

PLINK v1.90b6.24, GEMMA 0.98.5, GCTA 1.94.1, GxEMM 1.0 with LDAK 5, R 4.5.2 (cluster: R 4.5.3), Python 3.9.12.
R packages: rrBLUP 4.6.3, pROC 1.19.0.1, MASS, data.table, readxl, ggplot2, patchwork, dplyr, ggrepel, ggtext, ragg, scales.
