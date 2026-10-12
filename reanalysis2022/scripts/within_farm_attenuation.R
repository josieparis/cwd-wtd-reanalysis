# within-farm attenuation of the 2022 candidate SNPs (Results section 2): SNP-CWD correlation, marginal and within farm
# (dosage and CWD status both residualised on farm) for PRNP codon 96 and the non-PRNP candidates in Table S2
# retained = |r within farm| / |r marginal|, attenuation = 1 - retained

library(readxl)
library(tidyverse)

## working dir
setwd("<dir>/Seabury_data")

## read in covariates and candidate genotypes
cv <- read.csv("reanalysis2022/data/cov_1170.csv", stringsAsFactors = FALSE)
gt <- read.delim("reanalysis2022/data/candidate_genos_1170.tsv", check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(gt$CMSID == cv$CMSID))
y  <- suppressWarnings(as.numeric(cv$binary)); ok <- !is.na(y) & !is.na(cv$farm) & cv$farm != ""

## candidates from the authors' Table S2 (Dryad wh70rxwnt)
s2 <- suppressMessages(read_excel("doi_10_5061_dryad_wh70rxwnt__v20220622/TableS2_Revised_2022v2Final.xlsx"))
cand <- unique(sub("\\*+$", "", s2[["Marker*"]][grepl("^(Affx|G96S)", s2[["Marker*"]])]))   # strip footnote asterisks
miss <- setdiff(cand, names(gt)); cand <- intersect(cand, names(gt))

## marginal and within-farm correlation per candidate
att <- function(m) { g <- suppressWarnings(as.numeric(gt[[m]])); k <- ok & !is.na(g); f <- factor(cv$farm[k])
  rm <- cor(g[k], y[k]); rw <- cor(resid(lm(g[k] ~ f)), resid(lm(y[k] ~ f)))
  data.frame(marker = m, r_marginal = rm, r_within_farm = rw, retained = abs(rw) / abs(rm)) }
res <- map_dfr(cand, att)
pr <- res %>% filter(marker == "G96S_G_to_A_at286"); np <- res %>% filter(marker != "G96S_G_to_A_at286")

## summary
cat(sprintf("Table S2 candidates: %d unique markers (%d not among the deposited genotypes: %s)\n",
    length(cand) + length(miss), length(miss), paste(miss, collapse = ", ")))
cat(sprintf("PRNP: |r| %.3f -> %.3f within farm; attenuation %.1f%%\n",
    abs(pr$r_marginal), abs(pr$r_within_farm), 100 * (1 - pr$retained)))
cat(sprintf("non-PRNP candidates (n = %d): mean |r| %.3f -> %.3f; median attenuation %.0f%%; %d of %d lose > 50%%\n",
    nrow(np), mean(abs(np$r_marginal)), mean(abs(np$r_within_farm)), 100 * (1 - median(np$retained)),
    sum(np$retained < 0.5), nrow(np)))
