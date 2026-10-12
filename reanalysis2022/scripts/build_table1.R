# Table 1: lambda, number of SNPs at P < 5e-8 and P < 5e-5 and the PRNP P value
# for GEMMA, GCTA --mlma and GCTA --mlma-loco (2020 binary and score, 2022 binary)

## working dir
setwd("<dir>/Seabury_data")

PRNP <- "G96S_G_to_A_at286"

## genomic control lambda
lam  <- function(p) median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
## read p-values from an association output
rd   <- function(f, id, p) { d <- read.delim(f, stringsAsFactors = FALSE); setNames(as.numeric(d[[p]]), d[[id]]) }
## format P as e.g. 4.7 x 10^-11
fmtP <- function(p) { e <- floor(log10(p)); m <- p / 10^e
  if (round(m, 1) >= 10) { m <- m / 10; e <- e + 1 }
  sprintf("%.1f × 10^%d", m, e) }

## association results
## 2022 GEMMA uses its own centred GRM, 2022 GCTA uses the PLINK GRM (all 29 autosomes)
sets <- list(
  list(ds = "2020", n = 807, trait = "Binary", models = list(
    "GEMMA"            = rd("reanalysis/02_assoc/gemma_bin_results.tsv", "marker", "p_wald"),
    "GCTA --mlma"      = rd("reanalysis/02_assoc/output/gcta_mlma_bin.mlma", "SNP", "p"),
    "GCTA --mlma-loco" = rd("reanalysis/02_assoc/output/gcta_loco_bin.loco.mlma", "SNP", "p"))),
  list(ds = "2020", n = 807, trait = "Score", models = list(
    "GEMMA"            = rd("reanalysis/02_assoc/gemma_scr_results.tsv", "marker", "p_wald"),
    "GCTA --mlma"      = rd("reanalysis/02_assoc/output/gcta_mlma_scr.mlma", "SNP", "p"),
    "GCTA --mlma-loco" = rd("reanalysis/02_assoc/output/gcta__loco_scr.loco.mlma", "SNP", "p"))),
  list(ds = "2022", n = 1170, trait = "Binary", models = list(
    "GEMMA"            = rd("reanalysis2022/cluster_h2/gemma/output/pve_nocov.assoc.txt", "rs", "p_wald"),
    "GCTA --mlma"      = rd("reanalysis2022/mlma_results/mlma_bin_plink.mlma", "SNP", "p"),
    "GCTA --mlma-loco" = rd("reanalysis2022/mlma_results/loco_bin_chr29.loco.mlma", "SNP", "p"))))

## counts over the markers tested by all three models
rows <- list()
for (s in sets) {
  ok <- lapply(s$models, function(p) names(p)[is.finite(p) & p > 0])
  common <- Reduce(intersect, ok)
  cat(sprintf("%s %s: markers per model %s -> common %s\n", s$ds, s$trait,
              paste(sapply(ok, length), collapse = "/"), format(length(common), big.mark = ",")))
  for (m in names(s$models)) {
    p <- s$models[[m]][common]
    rows[[length(rows) + 1]] <- data.frame(
      Dataset = sprintf("%s (n = %s)", s$ds, format(s$n, big.mark = ",")), Phenotype = s$trait, Model = m,
      lambda = sprintf("%.2f", lam(p)), SNPs_P_5e8 = sum(p < 5e-8), SNPs_P_5e5 = sum(p < 5e-5),
      PRNP_P = fmtP(p[[PRNP]]), PRNP_P_raw = signif(p[[PRNP]], 4), markers = length(common),
      check.names = FALSE)
  }
}
t1 <- do.call(rbind, rows)
print(t1, row.names = FALSE)

## write table
dir.create("MANUSCRIPT_tables", showWarnings = FALSE)
write.csv(t1, "MANUSCRIPT_tables/table1_candidates.csv", row.names = FALSE)
cat("wrote MANUSCRIPT_tables/table1_candidates.csv\n")
