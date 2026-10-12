# Table 3: observed-scale SNP heritability (h2 +/- SE) of CWD in both datasets, for each covariate set and tool
# (original EMMAX value, GCTA --reml, GEMMA pve, rrBLUP); GCTA and rrBLUP use the PLINK --make-rel GRM
# also writes the long version read by fig4_heritability.R

library(tidyverse)

## working dir
setwd("<dir>/Seabury_data")

OUT <- "MANUSCRIPT_tables"
G20 <- "reanalysis2022/cluster_h2/gcta_2020/results"; G22 <- "reanalysis2022/cluster_h2/gcta_plink/results"
M20 <- "reanalysis/02_assoc/output";                  M22 <- "reanalysis2022/cluster_h2/gemma/output"

## h2, SE and n from a GCTA .hsq file
read_hsq <- function(f){ L <- strsplit(readLines(f), "\t")
  get <- function(key) as.numeric(Filter(function(z) z[1] == key, L)[[1]][-1])
  v <- get("V(G)/Vp"); c(h2 = v[1], se = v[2], n = get("n")[1]) }
## null-model pve, se(pve) and n from a GEMMA -lmm log (GEMMA builds its own GRM)
read_gemma <- function(f){ x <- readLines(f)
  num <- function(p) as.numeric(sub(".*= *", "", grep(p, x, value = TRUE)[1]))
  c(h2 = num("pve estimate in the null model"), se = num("se\\(pve\\) in the null model"),
    n = num("number of analyzed individuals")) }

## original 2020 values: the authors' deposited EMMAX log (Dryad xd2547dcw), one block per model
vc <- gsub("\r", "", readLines("doi_10_5061_dryad_xd2547dcw__v20200304/VarianceComponents807WTD.txt", warn = FALSE))
blocks <- split(vc, cumsum(grepl("^\\*{5,}", vc)))
emmax20 <- map_dfr(blocks, function(b){
  ph <- trimws(sub("Phenotype:", "", grep("^Phenotype:", b, value = TRUE)))
  if (!length(ph)) return(NULL)
  cv <- trimws(sub(" *\\(.*", "", grep("^  [A-Za-z_]+ \\([0-9]+\\)", b, value = TRUE)))
  val <- function(p) as.numeric(sub(".*: *", "", grep(p, b, value = TRUE)[1]))
  data.frame(pheno = ph, covs = paste(sort(cv), collapse = "+"),
             h2 = val("^Pseudo-heritability"), se = val("^Std Error of heritability"))
})
emx <- function(ph, cv = character(0)){
  r <- emmax20 %>% filter(pheno == ph, covs == paste(sort(cv), collapse = "+"))
  stopifnot(nrow(r) == 1); c(h2 = r$h2, se = r$se, n = 807) }
REG <- c("GeneralRegion", "Sex", "Age_OG")

## one row per model; original 2022 values as published (main text Fig 1a, Fig S1a,b)
rows <- list(
  list("2020","Binary","none",                              emx("CWD_Binary"),
       file.path(G20,"b_none.hsq"),         file.path(M20,"gwaa_bin_nocov.log.txt")),
  list("2020","Binary","+ sex/age/region",                  emx("CWD_Binary", REG),
       file.path(G20,"b_region.hsq"),       file.path(M20,"gwaa_bin_cov.log.txt")),
  list("2020","Binary","+ sex/age/region + tissue-score †", emx("CWD_Binary", c(REG, "CMS_CWDScore")),
       file.path(G20,"b_region_score.hsq"), file.path(M20,"gwaa_bin_covscore.log.txt")),
  list("2020","Score","none",                               emx("CMS_CWDScore"),
       file.path(G20,"s_none.hsq"),         file.path(M20,"gwaa_scr_nocov.log.txt")),
  list("2020","Score","+ sex/age/region",                   emx("CMS_CWDScore", REG),
       file.path(G20,"s_region.hsq"),       file.path(M20,"gwaa_scr_cov.log.txt")),
  list("2022","Binary","none",             c(h2 = 0.611, se = 0.056, n = 1170),
       file.path(G22,"h2_nocov.hsq"),       file.path(M22,"pve_nocov.log.txt")),
  list("2022","Binary","+ sex/age/region", c(h2 = 0.564, se = 0.059, n = 1166),
       file.path(G22,"h2_region.hsq"),      file.path(M22,"pve_region.log.txt")),
  list("2022","Binary","+ sex/age/farm",   c(h2 = 0.329, se = 0.067, n = 1151),
       file.path(G22,"h2_farm.hsq"),        file.path(M22,"pve_farm.log.txt")))

## rrBLUP estimates (h2_with_se.R) and check every tool used the same n
## 2020 n = 807; 2022 n = 1,170 (none), 1,166 (+ region), 1,151 (+ farm, farms with >= 10 deer)
rr <- read.csv("reanalysis2022/cluster_h2/rrblup/h2_rrblup_se.csv", stringsAsFactors = FALSE)
covkey <- function(x) sub(" †$", "", x)
long <- map_dfr(rows, function(r){
  b <- rr %>% filter(cohort == r[[1]], trait == tolower(r[[2]]), covset == covkey(r[[3]]))
  stopifnot(nrow(b) == 1)
  est <- rbind(`Original (EMMAX)` = r[[4]], GCTA = read_hsq(r[[5]]), GEMMA = read_gemma(r[[6]]),
               rrBLUP = c(h2 = b$h2, se = b$se, n = b$n))
  if (length(unique(est[, "n"])) != 1) stop("sample sizes differ in ", r[[1]], " ", r[[3]], ": ",
                                            paste(est[, "n"], collapse = "/"))
  data.frame(cohort = r[[1]], trait = r[[2]], covset = r[[3]], tool = rownames(est),
             n = est[, "n"], h2 = est[, "h2"], se = est[, "se"], row.names = NULL, check.names = FALSE)
})

## wide table
r3   <- function(x) sprintf("%.3f", round(x + 1e-9, 3))   # round half up
pm   <- function(h, s) paste(r3(h), "±", r3(s))
TOOL <- c("Original (EMMAX)", "GCTA", "GEMMA", "rrBLUP")
wide <- long %>%
  mutate(est = pm(h2, se), n = formatC(n, format = "d", big.mark = ",")) %>%
  pivot_wider(id_cols = c(cohort, trait, covset, n), names_from = tool, values_from = est) %>%
  select(Dataset = cohort, Trait = trait, Covariates = covset, n, all_of(TOOL))

## write tables
write.csv(wide, file.path(OUT, "table3_heritability_observed.csv"), row.names = FALSE)   # quoted: n has commas
write.csv(long, file.path(OUT, "table3_heritability_observed_long.csv"), row.names = FALSE)
print(as.data.frame(wide), row.names = FALSE, right = FALSE)
cat(sprintf("\nSE range: %.3f to %.3f\nwrote %s/table3_heritability_observed.csv (+ _long.csv)\n",
            min(long$se), max(long$se), OUT))
