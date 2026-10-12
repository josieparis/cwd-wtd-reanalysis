# genomic prediction of CWD status in 2022 (Table 2, Fig 2a,c): leave-one-farm-out (LOFO), a size-matched random
# control (20 replicates per farm) and random 5-fold CV, each without and with sex, age and region as covariates
# Rscript prediction_lofo.R [plink|custom] (default plink GRM), writes reanalysis2022/data/prediction_lofo_<grm>.rds

library(rrBLUP)
library(dplyr)

## working dir
setwd("<dir>/Seabury_data")

GRM <- commandArgs(TRUE)[1]; if (is.na(GRM)) GRM <- "plink"

## read in covariates and GRM
cov  <- read.csv("reanalysis2022/data/cov_1170.csv", stringsAsFactors = FALSE)
n    <- nrow(cov); ids <- cov$CMSID
if (GRM == "plink") {
  K <- as.matrix(read.table(gzfile("reanalysis2022/01_structure/grm2022.rel.gz")))
  kid <- read.table("reanalysis2022/01_structure/grm2022.rel.id", stringsAsFactors = FALSE)$V2
  dimnames(K) <- list(kid, kid); K <- K[ids, ids]
} else {
  K <- matrix(readBin("reanalysis2022/data/K_1170.bin", "double", n * n, size = 8), n, n)   # old in-house GRM
}
## bend to positive definite
ei <- eigen(K, symmetric = TRUE); K <- ei$vectors %*% (pmax(ei$values, 1e-6) * t(ei$vectors))
dimnames(K) <- list(ids, ids)

## phenotype and covariates (none, or sex + age + region); row = row of the deer in K
cov  <- cov %>% mutate(row = row_number(), y = suppressWarnings(as.numeric(binary)),
                       sexn = suppressWarnings(as.numeric(sex)))
y    <- cov$y; farm <- cov$farm
cc   <- cov %>% filter(!is.na(sexn), !is.na(age))                   # deer with known sex and age
Xall <- list(none = matrix(1, n, 1),
             cov  = { X <- matrix(NA_real_, n, 5)
                      X[cc$row, ] <- model.matrix(~ factor(sexn) + age + region, cc); X })

## GBLUP on the training animals, prediction = fixed effects + genomic value
gblup <- function(train, X) {                       # predicted phenotype for every animal
  tr <- train[complete.cases(X[train, , drop = FALSE])]
  Z  <- matrix(0, length(tr), n); Z[cbind(seq_along(tr), tr)] <- 1
  fit <- mixed.solve(y = y[tr], Z = Z, K = K, X = X[tr, , drop = FALSE])
  as.numeric(X %*% fit$beta) + as.numeric(fit$u)
}
## AUC, accuracy, sensitivity, specificity, MCC, RMSE and r, with cases called at >= 0.5 (authors' rule)
metvec <- function(p, t, thr = 0.5) {
  ok <- !is.na(p) & !is.na(t); p <- p[ok]; t <- t[ok]
  np <- sum(t == 1); nn <- sum(t == 0); N <- length(t)
  if (!np || !nn) return(setNames(rep(NA, 7), MET))
  auc <- (sum(rank(p)[t == 1]) - np * (np + 1) / 2) / (np * nn)
  cls <- as.integer(p >= thr)
  TP <- sum(cls == 1 & t == 1); TN <- sum(cls == 0 & t == 0); FP <- sum(cls == 1 & t == 0); FN <- sum(cls == 0 & t == 1)
  den <- sqrt(as.double(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  c(auc = auc, acc = (TP + TN) / N, sens = TP / (TP + FN), spec = TN / (TN + FP),
    mcc = if (den > 0) (TP * TN - FP * FN) / den else 0, rmse = sqrt(mean((p - t)^2)),
    r = if (sd(p) > 0) cor(p, t) else NA_real_)
}
MET <- c("auc", "acc", "sens", "spec", "mcc", "rmse", "r")

## deer with known CWD status on a named farm, and farms with >= 10 deer and >= 3 cases and controls
vd     <- cov %>% filter(!is.na(y), !is.na(farm), farm != "")
valid  <- vd$row
qfarms <- vd %>% group_by(farm) %>% summarise(n = n(), cases = sum(y == 1), controls = sum(y == 0), .groups = "drop") %>%
  filter(n >= 10, cases >= 3, controls >= 3) %>% pull(farm)
REPS   <- 20

## run LOFO, size-matched and random 5-fold (seeds as in the original leave_one_farm_out.R)
out <- list()
for (cs in names(Xall)) {
  X <- Xall[[cs]]; v <- valid[complete.cases(X[valid, , drop = FALSE])]
  lofo <- size_m <- matrix(NA, length(qfarms), 7, dimnames = list(qfarms, MET))
  for (fi in seq_along(qfarms)) {
    f <- qfarms[fi]; fidx <- v[farm[v] == f]; m <- length(fidx)
    lofo[fi, ] <- metvec(gblup(setdiff(v, fidx), X)[fidx], y[fidx])
    rm <- t(sapply(1:REPS, function(r) { set.seed(1000L * fi + r); ti <- sample(v, m)
      metvec(gblup(setdiff(v, ti), X)[ti], y[ti]) }))
    size_m[fi, ] <- colMeans(rm, na.rm = TRUE)
    cat(sprintf("[%s] %-9s n=%3d | LOFO auc=%.3f sens=%.3f mcc=%+.3f | size-matched auc=%.3f\n",
                cs, f, m, lofo[fi, "auc"], lofo[fi, "sens"], lofo[fi, "mcc"], size_m[fi, "auc"]))
  }
  set.seed(42); fold <- sample(rep(1:5, length.out = length(v)))
  cv <- t(sapply(1:5, function(k) { ti <- v[fold == k]; metvec(gblup(setdiff(v, ti), X)[ti], y[ti]) }))
  out[[cs]] <- list(lofo = lofo, size_matched = size_m, random_5fold = cv, n = length(v))
}
out$grm <- GRM; out$farms <- qfarms

## save
f <- sprintf("reanalysis2022/data/prediction_lofo_%s.rds", GRM)
saveRDS(out, f)

## summary, mean +/- SD
ms <- function(m) sprintf("%.3f ± %.3f", colMeans(m, na.rm = TRUE), apply(m, 2, sd, na.rm = TRUE))
for (cs in c("none", "cov")) {
  cat(sprintf("\n=== %s (GRM: %s; n = %d) ===\n", cs, GRM, out[[cs]]$n))
  print(data.frame(metric = MET, random_5fold = ms(out[[cs]]$random_5fold),
                   size_matched = ms(out[[cs]]$size_matched), LOFO = ms(out[[cs]]$lofo)), row.names = FALSE)
  cat("LOFO farms with sensitivity = 0:", sum(out[[cs]]$lofo[, "sens"] == 0, na.rm = TRUE), "of", length(qfarms), "\n")
}
cat("\nwrote", f, "\n")
