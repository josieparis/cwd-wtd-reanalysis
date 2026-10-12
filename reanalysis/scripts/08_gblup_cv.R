# Step 3a: reproduce Seabury's GBLUP prediction accuracy with rrBLUP (Endelman 2011) and our GRM
# random 5-fold CV x 50 iterations; binary = AUC of the pooled held-out predictions, score = Pearson r
# compare with AF6/AF8 and AF10/AF12

library(rrBLUP)
library(pROC)

## working dir
ROOT <- "<dir>/Seabury_data/reanalysis"
S <- file.path(ROOT,"01_structure")

## read in GRM
ids <- read.table(file.path(S,"grm.rel.id"))[,1]
n <- length(ids)
K <- as.matrix(read.table(gzfile(file.path(S,"grm.rel.gz"))))
dimnames(K) <- list(ids, ids)
## nudge to positive definite for the solver
K <- K + diag(1e-4, n)

## phenotypes
meta <- read.delim(file.path(ROOT,"00_data/samples.tsv"), colClasses=c(clean_id="character"))
rownames(meta) <- meta$clean_id; meta <- meta[ids,]
ybin <- meta$binary; yscr <- meta$score

set.seed(20200304)
Kfold <- 5; ITER <- 50

## one round of k-fold CV
cv_once <- function(y, is_binary){
  fold <- sample(rep(1:Kfold, length.out=n))
  pred <- rep(NA_real_, n)
  for(k in 1:Kfold){
    test <- which(fold==k)
    df <- data.frame(gid=ids, yy=y)
    df$yy[test] <- NA
    fit <- kin.blup(data=df, geno="gid", pheno="yy", K=K)
    pred[test] <- fit$g[test]
  }
  if(is_binary) as.numeric(pROC::auc(response=y, predictor=pred, quiet=TRUE, levels=c(0,1), direction="<"))
  else cor(pred, y)
}

## repeat CV ITER times
run <- function(y, is_binary, label){
  v <- numeric(ITER)
  for(i in 1:ITER){ v[i] <- cv_once(y, is_binary); if(i%%10==0) cat(sprintf("  %s iter %d\n",label,i)) }
  cat(sprintf(">>> %s: mean=%.3f  sd=%.3f  (min %.3f max %.3f)\n", label, mean(v), sd(v), min(v), max(v)))
  v
}

## run
cat("=== GBLUP 5-fold CV x50, our GRM ===\n")
ab <- run(ybin, TRUE,  "BINARY AUC")
as <- run(yscr, FALSE, "SCORE  r  ")

## save
write.csv(data.frame(iter=1:ITER, binary_auc=ab, score_r=as),
          file.path(S,"gblup_cv_random.csv"), row.names=FALSE)
cat("\nReference (Seabury AF): binary AUC~0.831-0.833, score r~0.584-0.588\n")
cat("saved gblup_cv_random.csv\n")
