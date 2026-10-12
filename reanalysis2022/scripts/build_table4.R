# Table 4: liability-scale heritability of CWD at K = 0.05, 0.10 and 0.26, plus its footnotes
# GCTA-GREML + Lee (GCTA --prevalence, no covariates) and LDAK-PCGC (GxEMM HOM model, binary = TRUE),
# both datasets, all on the PLINK GRM

## working dir
setwd("<dir>/Seabury_data")

OUT <- "MANUSCRIPT_tables"
KS  <- c("0.05", "0.10", "0.26")
## Seabury et al. (2022) only give the K = 0.05 value in the text (they fitted K = 0.01-0.35)
SEABURY <- c("0.05" = "0.557 ± 0.053", "0.10" = "n.r.", "0.26" = "n.r.")

pm   <- function(h, s) sprintf("%.3f ± %.3f", h, s)
cells <- function(h, s) {                                       # dagger if > 1, section sign if the SE collapses
  out <- paste0(pm(h, s), ifelse(h > 1, "†", ""))
  if (s[2] < 0.5 * min(s[1], s[3])) out[2] <- paste0(out[2], "§")
  out
}
read_hsq <- function(f) {
  L <- strsplit(readLines(f), "\t")
  get <- function(key) as.numeric(Filter(function(z) z[1] == key, L)[[1]][-1])
  list(obs = get("V(G)/Vp"), liab = get("V(G)/Vp_L"), n = get("n")[1])
}
hom <- function(o) c(h = unname(o$h2[1]), s = sqrt(o$h2Covmat[1, 1]))

## GCTA-GREML + Lee (no covariates)
g22 <- lapply(KS, function(k) read_hsq(sprintf("reanalysis2022/cluster_h2/gcta_plink/results/h2_nocov_liab_%s.hsq", k)))
g20 <- lapply(KS, function(k) read_hsq(sprintf("reanalysis2022/cluster_h2/gcta_2020/results/b_none_liab_%s.hsq", k)))
gl  <- function(g) list(h = sapply(g, function(x) x$liab[1]), s = sapply(g, function(x) x$liab[2]))

## PCGC (GxEMM HOM, binary = TRUE)
## 2022 region and farm rows from the farm-adjusted run (same 1,120 deer), 2020 from the 2020 run
fa  <- readRDS("reanalysis2022/cluster_gxemm/gxemm_farmadj_results.rds")
t20 <- readRDS("reanalysis2022/cluster_gxemm_2020/gxemm_2020_results.rds")
pg  <- function(objs) { m <- sapply(objs, hom); list(h = m["h", ], s = m["s", ]) }
pA  <- pg(lapply(KS, function(k) fa$A_region_adj[[paste0("liab_", k)]]$hom))
pB  <- pg(lapply(KS, function(k) fa$B_farm_adj[[paste0("liab_", k)]]$hom))
p20 <- pg(lapply(KS, function(k) t20[[paste0("liab_hom_", k)]]))
c807 <- read.csv("reanalysis/01_structure/cov_807.csv", stringsAsFactors = FALSE)
num  <- function(x) suppressWarnings(as.numeric(x))
n20p <- sum(!is.na(num(c807$binary)) & !is.na(num(c807$sex)) & !is.na(num(c807$age)))   # same filter as run_gxemm_2020.R

## build and write table
r22 <- gl(g22); r20 <- gl(g20)
nf  <- function(n) formatC(as.integer(n), format = "d", big.mark = ",")
rows <- rbind(
  c("GCTA-GREML + Lee", "2022", "Seabury et al. 2022", "not stated", nf(1170), SEABURY),
  c("GCTA-GREML + Lee", "2022", "this study", "none", nf(g22[[1]]$n), cells(r22$h, r22$s)),
  c("GCTA-GREML + Lee", "2020", "this study", "none", nf(g20[[1]]$n), cells(r20$h, r20$s)),
  c("LDAK-PCGC", "2022", "this study", "region + sex + age", nf(fa$meta$N), cells(pA$h, pA$s)),
  c("LDAK-PCGC", "2022", "this study", "farm + sex + age",   nf(fa$meta$N), cells(pB$h, pB$s)),
  c("LDAK-PCGC", "2020", "this study", "region + sex + age", nf(n20p),      cells(p20$h, p20$s))
)
colnames(rows) <- c("Estimator", "Dataset", "Source", "Covariates", "n", paste("K =", KS))
write.csv(rows, file.path(OUT, "table4_heritability_liability.csv"), row.names = FALSE)   # quoted: n has commas

## footnotes
o22 <- pm(g22[[1]]$obs[1], g22[[1]]$obs[2]); o20 <- pm(g20[[1]]$obs[1], g20[[1]]$obs[2])
notes <- c(
  "**Table 4.** Liability-scale heritability of CWD status at assumed population prevalence K (h² ± SE).",
  "K = 0.26 is the weighted-mean farm prevalence reported by Seabury et al. (2022). The rows are two",
  "different estimators and are not like-for-like: compare within an estimator (e.g. region- vs",
  "farm-adjusted PCGC), not across.", "",
  paste0("- **GCTA-GREML + Lee:** observed-scale REML heritability transformed with the Lee et al. (2011) ",
         "conversion (GCTA `--prevalence`), the method used by Seabury et al. (2022). Our rows transform the ",
         "no-covariate observed estimates (2022: ", o22, "; 2020: ", o20, "; Table 3)."),
  paste0("- **PCGC:** phenotype-correlation genotype-correlation regression (Golan et al. 2014), fit as the ",
         "homogeneous (HOM) GxEMM model with `binary = TRUE` (Dahl et al. 2020), i.e. a single-GRM PCGC run in ",
         "LDAK after adjusting the GRM for the covariates, on the PLINK 1.9 `--make-rel` VanRaden GRM; SEs are ",
         "jackknife-based."),
  paste0("- 2022 PCGC rows use the ", formatC(fa$meta$N, big.mark = ","), " deer on the ", length(fa$meta$farms), " farms with ≥ 10 deer ",
         "and both cases and controls (complete sex and age), so the region- and farm-adjusted rows are on ",
         "identical animals."),
  paste0("- n.r., not reported: Seabury et al. (2022) fitted K = 0.01–0.35 but give only the K = 0.05 ",
         "value in the text."),
  "- † Inadmissible (heritability > 1).",
  "- § SE anomalously small relative to the adjacent prevalences, where the estimate crosses 1; interpret with caution."
)
writeLines(notes, file.path(OUT, "table4_heritability_liability_notes.md"))
cat("wrote", file.path(OUT, "table4_heritability_liability.csv"), "+ notes\n")
