# Table 2: genomic prediction accuracy in the 2022 dataset (PLINK GRM, from prediction_lofo.R)
# mean +/- SD over the 5 folds, or over the 13 held-out farms (size-matched: mean of 20 replicates per farm)
# extra row: accuracy of calling every deer CWD-negative on the same test sets
# also prints the numbers quoted in the text

## working dir
setwd("<dir>/Seabury_data")

## read in prediction results
P   <- readRDS("reanalysis2022/data/prediction_lofo_plink.rds")
MET <- c(auc = "AUC", acc = "Accuracy", neg = "Accuracy if every deer is called CWD-negative",
         sens = "Sensitivity", spec = "Specificity", mcc = "MCC", rmse = "RMSE", r = "Predictive r")
DES <- c(random_5fold = "Random 5-fold", size_matched = "Sample size-matched", lofo = "Leave-one-farm-out")

## accuracy of calling every deer CWD-negative (= share of controls)
## same test sets as prediction_lofo.R (same animals and seeds)
cv    <- read.csv("reanalysis2022/data/cov_1170.csv", stringsAsFactors = FALSE)
y     <- suppressWarnings(as.numeric(cv$binary)); farm <- cv$farm; sexn <- suppressWarnings(as.numeric(cv$sex))
valid <- which(!is.na(y) & !is.na(farm) & farm != "")
neg   <- function(t) if (sum(t == 1) && sum(t == 0)) mean(t == 0) else NA
for (cs in c("none", "cov")) {
  v <- if (cs == "none") valid else valid[!is.na(sexn[valid]) & !is.na(cv$age[valid])]
  stopifnot(length(v) == P[[cs]]$n)
  P[[cs]]$lofo <- cbind(P[[cs]]$lofo, neg = sapply(P$farms, function(f) neg(y[v[farm[v] == f]])))
  P[[cs]]$size_matched <- cbind(P[[cs]]$size_matched, neg = sapply(seq_along(P$farms), function(fi) {
    m <- sum(farm[v] == P$farms[fi])
    mean(sapply(1:20, function(r) { set.seed(1000L * fi + r); neg(y[sample(v, m)]) }), na.rm = TRUE) }))
  set.seed(42); fold <- sample(rep(1:5, length.out = length(v)))
  P[[cs]]$random_5fold <- cbind(P[[cs]]$random_5fold, neg = sapply(1:5, function(k) neg(y[v[fold == k]])))
}

## mean +/- SD (2 dp) with a proper minus sign
ms  <- function(m, k) { x <- m[, k]; s <- sprintf("%.2f ± %.2f", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)); sub("^-", "−", s) }

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

## accuracy vs calling every deer CWD-negative (no covariates)
aL <- P$none$lofo[, "acc"]; bL <- P$none$lofo[, "neg"]
cat(sprintf("\nAccuracy vs calling every deer negative: random 5-fold %.2f vs %.2f | leave-one-farm-out %.2f vs %.2f (model below the rule on %d of %d farms, tied on %d)\n",
    mean(P$none$random_5fold[, "acc"]), mean(P$none$random_5fold[, "neg"]), mean(aL), mean(bL),
    sum(aL < bL - 1e-9), length(aL), sum(abs(aL - bL) < 1e-9)))
cat("\nwrote MANUSCRIPT_tables/table2_prediction.csv\n")
