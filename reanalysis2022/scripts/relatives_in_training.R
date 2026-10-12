# how many test deer have a relative in the training set under leave-one-farm-out (the 13 farms of Table 2)
# vs random 5-fold CV (10 replicates, as in relatedness_binned_auc.R), from the max PLINK GRM value (Results section 3)

## working dir
setwd("<dir>/Seabury_data")

source("reanalysis2022/scripts/load_grm.R")

## read in covariates and PLINK GRM
cv <- read.csv("reanalysis2022/data/cov_1170.csv", stringsAsFactors = FALSE)
y  <- suppressWarnings(as.numeric(cv$binary))
K  <- load_grm2022(cv$CMSID, bend = FALSE); diag(K) <- NA

## print median max GRM and % of test deer with a first-degree (>= 0.4) or close (>= 0.2) relative in training
summ <- function(lab, mx) cat(sprintf("%-22s n = %5d | median max GRM %.3f | first-degree (>= 0.4) %.1f%% | close (>= 0.2) %.1f%%\n",
                                      lab, length(mx), median(mx), 100 * mean(mx >= 0.4), 100 * mean(mx >= 0.2)))

## leave-one-farm-out (farms with >= 10 deer and >= 3 cases and controls)
valid <- which(!is.na(y) & !is.na(cv$farm) & cv$farm != "")
tab <- table(cv$farm[valid], y[valid]); q <- rownames(tab)[rowSums(tab) >= 10 & apply(tab, 1, min) >= 3]
summ("leave-one-farm-out", unlist(lapply(q, function(f) { te <- valid[cv$farm[valid] == f]
  apply(K[te, setdiff(valid, te), drop = FALSE], 1, max) })))

## random 5-fold, 10 replicates
ph <- which(!is.na(y)); set.seed(1)                      # seed as in relatedness_binned_auc.R
summ("random 5-fold (x10)", unlist(lapply(1:10, function(it) { fold <- sample(rep(1:5, length.out = length(ph)))
  unlist(lapply(1:5, function(k) apply(K[ph[fold == k], ph[fold != k], drop = FALSE], 1, max))) })))
