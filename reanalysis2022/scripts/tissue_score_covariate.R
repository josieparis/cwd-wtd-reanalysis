# the 2020 model that adjusts binary CWD for the ordinal tissue score
# (Results section 2, Discussion and Supplementary Fig 2d)

library(tidyverse)

## working dir
setwd("<dir>/Seabury_data")

## 1. binary = score > 0 (every control scores 0, every case 1 or 2)
cv <- read.csv("reanalysis/01_structure/cov_807.csv", stringsAsFactors = FALSE)
y <- cv$binary; s <- cv$score
stopifnot(all(s[y == 0] == 0), all(s[y == 1] >= 1))
cat(sprintf("1. controls %d (all score 0) | cases %d (score 1: %d, score 2: %d)\n",
            sum(y == 0), sum(y == 1), sum(s == 1), sum(s == 2)))

## 2. variance explained by the fixed covariates in the authors' EMMAX log (Dryad xd2547dcw), without and with the score
vc <- gsub("\r", "", readLines("doi_10_5061_dryad_xd2547dcw__v20200304/VarianceComponents807WTD.txt", warn = FALSE))
blk <- split(vc, cumsum(grepl("^\\*{5,}", vc)))
fx <- map_dfr(blk, function(b) {
  ph <- trimws(sub("Phenotype:", "", grep("^Phenotype:", b, value = TRUE)))
  pf <- grep("^Prop. explained by fixed covariates", b, value = TRUE)
  if (!length(ph) || !length(pf)) return(NULL)
  data.frame(pheno = ph, score_covariate = any(grepl("CMS_CWDScore \\(", b)),
             prop_fixed = as.numeric(sub(".*: *", "", pf)))
})
cat("2. authors' EMMAX log, proportion of phenotypic variance explained by the fixed covariates:\n")
print(fx, row.names = FALSE)

## 3. residual of binary ~ region + sex + age + score, regressed on lymph node only (score == 1)
r  <- resid(lm(binary ~ factor(region) + factor(sex) + age + score, data = cv))
r2 <- summary(lm(r ~ I(s == 1)))$r.squared
cat(sprintf("3. residual of binary ~ region + sex + age + score explained by the lymph-node-only contrast: %.1f%%\n", 100 * r2))

## 4. GEMMA single-GRM scans (autosomal SNPs, as in Table 1) with and without the score, and on the score itself,
## for PRNP and the two strongest non-PRNP markers in the authors' Table S1
auto <- read.table("reanalysis/02_assoc/deer_bin.bim", stringsAsFactors = FALSE)$V2
gm <- function(f) { d <- read.delim(file.path("reanalysis/02_assoc/output", f), stringsAsFactors = FALSE) %>%
  filter(rs %in% auto); setNames(d$p_wald, d$rs) }
bc <- gm("gwaa_bin_cov.assoc.txt"); bs <- gm("gwaa_bin_covscore.assoc.txt"); sc <- gm("gwaa_scr_cov.assoc.txt")
S <- c(PRNP = "G96S_G_to_A_at286", `NOVA1-FOXG1 (chr21)` = "Affx-573993815", `near KLF6 (chr13)` = "Affx-574078732")
cat("4. GEMMA, binary ~ sex + age + region [+ score]; and score ~ sex + age + region\n")
for (k in names(S)) { m <- S[[k]]
  cat(sprintf("   %-20s binary: P = %.2g (rank %d of %d) -> with score P = %.2g (rank %d) | tissue score: P = %.2g\n",
      k, bc[[m]], rank(bc)[[m]], length(bc), bs[[m]], rank(bs)[[m]], sc[[m]])) }
cat(sprintf("   genome-wide significant (P < 5e-8) with the score: %s\n", paste(names(bs)[bs < 5e-8], collapse = ", ")))
