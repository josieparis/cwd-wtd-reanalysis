# Step 3b: is the GBLUP AUC of 0.85 real genetics or structure? same GBLUP and GRM as 3a
# A) random 5-fold, B) AUC by each test animal's max relatedness to training, C) leave-one-region-out,
# D) PRNP G96S only, E) region only, F) PRNP + region

library(rrBLUP)
library(pROC)

## working dir
ROOT <- "<dir>/Seabury_data/reanalysis"
S <- file.path(ROOT,"01_structure")
auc <- function(y,p) as.numeric(pROC::auc(response=y, predictor=p, quiet=TRUE, levels=c(0,1), direction="<"))

## read in GRM and phenotypes
ids <- read.table(file.path(S,"grm.rel.id"))[,1]; n <- length(ids)
K <- as.matrix(read.table(gzfile(file.path(S,"grm.rel.gz")))); dimnames(K) <- list(ids,ids)
Kd <- K + diag(1e-4,n)
meta <- read.delim(file.path(ROOT,"00_data/samples.tsv"), colClasses=c(clean_id="character"))
rownames(meta) <- meta$clean_id; meta <- meta[ids,]
y <- meta$binary; region <- meta$region
## PRNP genotypes (plink --recode A, see 00b_plink_2020.sh)
g96 <- read.table(file.path(ROOT, "work", "prnp.raw"), header=TRUE)
rownames(g96) <- g96$IID; g96 <- g96[ids,]; G96 <- g96$G96S_G_to_A_at286_A

## GBLUP prediction of the test animals
gblup_pred <- function(train_idx, test_idx, yy){
  df <- data.frame(gid=ids, yv=yy); df$yv[test_idx] <- NA
  kin.blup(data=df, geno="gid", pheno="yv", K=Kd)$g[test_idx]
}

set.seed(20200304); ITER <- 30; Kfold <- 5

## A + B: random CV, keep predictions and each test animal's max relatedness to training
pred_all <- c(); y_all <- c(); maxrel_all <- c()
Koff <- K; diag(Koff) <- -Inf
for(it in 1:ITER){
  fold <- sample(rep(1:Kfold,length.out=n))
  for(k in 1:Kfold){
    te <- which(fold==k); tr <- which(fold!=k)
    p <- gblup_pred(tr,te,y)
    mr <- apply(Koff[te,tr,drop=FALSE],1,max)   # max GRM to any training animal
    pred_all<-c(pred_all,p); y_all<-c(y_all,y[te]); maxrel_all<-c(maxrel_all,mr)
  }
}
cat(sprintf("A) random 5-fold GBLUP AUC (pooled, %d iter) = %.3f\n", ITER, auc(y_all,pred_all)))
cat("\nB) AUC stratified by test animal's MAX relatedness to training set:\n")
brk <- c(-Inf,0.05,0.125,0.25,Inf); lab <- c("<0.05 (unrelated)","0.05-0.125 (~3rd-4th deg)","0.125-0.25 (~2nd deg)",">0.25 (1st-2nd deg)")
bin <- cut(maxrel_all,brk,labels=lab)
for(l in lab){ m<-bin==l; if(sum(m)>20) cat(sprintf("   %-26s n=%5d  AUC=%.3f  prevalence=%.2f\n",l,sum(m),auc(y_all[m],pred_all[m]),mean(y_all[m]))) }

## C: leave-one-region-out GBLUP
cat("\nC) leave-one-region-out GBLUP (train other 2 regions, predict held-out region):\n")
for(r in c("South","Midwest","East")){
  te <- which(region==r); tr <- which(region!=r)
  p <- gblup_pred(tr,te,y)
  cat(sprintf("   predict %-8s n=%3d prev=%.2f  AUC=%.3f\n", r, length(te), mean(y[te]), auc(y[te],p)))
}

## D: PRNP G96S only (5-fold logistic, 30 iterations)
av <- numeric(ITER)
for(it in 1:ITER){
  fold <- sample(rep(1:Kfold,length.out=n)); pr <- rep(NA,n)
  for(k in 1:Kfold){ te<-which(fold==k); tr<-which(fold!=k)
    m<-glm(y[tr]~G96[tr],family=binomial); pr[te]<-predict(m,data.frame(G96=G96[te]),type="response") }
  av[it]<-auc(y,pr)
}
cat(sprintf("\nD) PRNP G96S-only predictor AUC = %.3f (sd %.3f)\n", mean(av), sd(av)))

## E: region only (no genotypes)
ar <- numeric(ITER)
for(it in 1:ITER){
  fold <- sample(rep(1:Kfold,length.out=n)); pr <- rep(NA,n)
  for(k in 1:Kfold){ te<-which(fold==k); tr<-which(fold!=k)
    tab<-tapply(y[tr],region[tr],mean); pr[te]<-tab[region[te]] }
  ar[it]<-auc(y,pr)
}
cat(sprintf("E) region-only predictor AUC = %.3f (sd %.3f)\n", mean(ar), sd(ar)))

## F: PRNP + region (no genome-wide GRM)
af <- numeric(ITER)
for(it in 1:ITER){
  fold <- sample(rep(1:Kfold,length.out=n)); pr <- rep(NA,n)
  for(k in 1:Kfold){ te<-which(fold==k); tr<-which(fold!=k)
    m<-glm(y[tr]~G96[tr]+factor(region[tr]),family=binomial)
    pr[te]<-predict(m,data.frame(`G96[tr]`=G96[te],`region[tr]`=region[te],check.names=FALSE),type="response") }
  af[it]<-auc(y,pr)
}
cat(sprintf("F) PRNP + region predictor AUC = %.3f (sd %.3f)\n", mean(af), sd(af)))
cat("\nReference: full-GRM random-CV GBLUP AUC ~0.85 (3a)\n")
