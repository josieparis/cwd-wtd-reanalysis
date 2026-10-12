# Table 2: genomic prediction accuracy in the 2022 dataset (PLINK GRM, from prediction_lofo.R)
# mean +/- SD over the 5 folds, or over the 13 held-out farms (size-matched: mean of 20 replicates per farm)
# also prints the numbers quoted in the text

## working dir
setwd("<dir>/Seabury_data")

## read in prediction results
P   <- readRDS("reanalysis2022/data/prediction_lofo_plink.rds")
MET <- c(auc = "AUC", acc = "Accuracy", sens = "Sensitivity", spec = "Specificity", mcc = "MCC",
         rmse = "RMSE", r = "Predictive r")
DES <- c(random_5fold = "Random 5-fold", size_matched = "Sample size-matched", lofo = "Leave-one-farm-out")
## mean +/- SD with a proper minus sign
ms  <- function(m, k) { x <- m[, k]; s <- sprintf("%.3f ± %.3f", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)); sub("^-", "−", s) }

## build table
tab <- data.frame(Metric = unname(MET), check.names = FALSE)
for (cs in c("none", "cov")) for (d in names(DES))
  tab[[sprintf("%s | %s", c(none = "No covariates", cov = "Covariates (+ sex / age / region)")[[cs]], DES[[d]])]] <-
    sapply(names(MET), function(k) ms(P[[cs]][[d]], k))

## write table
dir.create("MANUSCRIPT_tables", showWarnings = FALSE)
write.csv(tab, "MANUSCRIPT_tables/table2_prediction.csv", row.names = FALSE)
print(tab, row.names = FALSE)

## LOFO farms with zero sensitivity and the false-negative rate
for (cs in c("none", "cov")) {
  L <- P[[cs]]$lofo
  cat(sprintf("\n[%s] n = %d deer, %d farms | LOFO farms with sensitivity = 0: %d (%s) | mean false-negative rate %.2f\n",
      cs, P[[cs]]$n, nrow(L), sum(L[, "sens"] == 0, na.rm = TRUE),
      paste(rownames(L)[which(L[, "sens"] == 0)], collapse = ", "), 1 - mean(L[, "sens"], na.rm = TRUE)))
}

## accuracy vs calling every deer CWD-negative (accuracy = share of controls)
## same test sets as Table 2 (no covariates, folds as in prediction_lofo.R, seed 42)
cv <- read.csv("reanalysis2022/data/cov_1170.csv", stringsAsFactors = FALSE)
y  <- suppressWarnings(as.numeric(cv$binary)); farm <- cv$farm
valid <- which(!is.na(y) & !is.na(farm) & farm != ""); stopifnot(length(valid) == P$none$n)
set.seed(42); fold <- sample(rep(1:5, length.out = length(valid)))
b5 <- sapply(1:5, function(k) mean(y[valid[fold == k]] == 0))
bL <- sapply(rownames(P$none$lofo), function(f) mean(y[valid[farm[valid] == f]] == 0))
aL <- P$none$lofo[, "acc"]
cat(sprintf("\nAccuracy vs calling every deer negative: random 5-fold %.3f vs %.3f | leave-one-farm-out %.3f vs %.3f (model below the rule on %d of %d farms, tied on %d)\n",
    mean(P$none$random_5fold[, "acc"]), mean(b5), mean(aL), mean(bL),
    sum(aL < bL - 1e-9), length(aL), sum(abs(aL - bL) < 1e-9)))
cat("\nwrote MANUSCRIPT_tables/table2_prediction.csv\n")
