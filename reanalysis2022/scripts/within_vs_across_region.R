# within-region vs across-region genomic prediction (Fig 2d): 5-fold CV inside each region, leave-one-region-out,
# and leave-one-region-out with the top 5 genomic PCs as covariates (design after Fritsche-Neto, Mol Breeding 2018)
# writes reanalysis2022/data/within_vs_across_region_summary.csv and within_vs_across_region_perregion.csv

library(rrBLUP)

## working dir
setwd("<dir>/Seabury_data")

source("reanalysis2022/scripts/load_grm.R")

## read in data
cov <- read.csv("reanalysis2022/data/cov_1170.csv", stringsAsFactors=FALSE)
n   <- nrow(cov)
y   <- suppressWarnings(as.numeric(cov$binary))
reg <- cov$region
ids <- cov$CMSID

## PLINK GRM and top 5 genomic PCs
Kp <- load_grm2022(ids)                                     # PLINK GRM, bent to positive definite
ei <- eigen(Kp, symmetric=TRUE)
PC <- ei$vectors[,1:5]; colnames(PC) <- paste0("PC",1:5)   # for the PC covariate test

## GBLUP prediction: intercept + genomic value (+ PC fixed effects)
gblup <- function(ymask, usePC=FALSE){
  d <- data.frame(gid=ids, yy=ymask, PC)
  if(usePC) fit <- kin.blup(d, geno="gid", pheno="yy", K=Kp, covariate=paste0("PC",1:5))
  else      fit <- kin.blup(d, geno="gid", pheno="yy", K=Kp)
  g <- fit$g[match(ids, names(fit$g))]
  mean(ymask, na.rm=TRUE) + g
}

## metrics as in Seabury 2020 Table 1, with cases called at >= 0.5
metvec <- function(pred, truth, thr=0.5){
  ok <- !is.na(pred) & !is.na(truth); p <- pred[ok]; t <- truth[ok]
  npos <- sum(t==1); nneg <- sum(t==0); N <- length(t)
  if(!npos || !nneg) return(c(auc=NA,acc=NA,sens=NA,spec=NA,mcc=NA,rmse=NA,r=NA))
  rk <- rank(p); auc <- (sum(rk[t==1]) - npos*(npos+1)/2)/(npos*nneg)
  rmse <- sqrt(mean((p - t)^2)); r <- if(sd(p)>0) cor(p,t) else NA_real_
  cls <- as.integer(p >= thr)
  TP<-sum(cls==1&t==1); TN<-sum(cls==0&t==0); FP<-sum(cls==1&t==0); FN<-sum(cls==0&t==1)
  acc<-(TP+TN)/N; sens<-TP/(TP+FN); spec<-TN/(TN+FP)
  den<-sqrt(as.double(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)); mcc<-if(den>0)(TP*TN-FP*FN)/den else 0
  c(auc=auc,acc=acc,sens=sens,spec=spec,mcc=mcc,rmse=rmse,r=r)
}
MET <- c("auc","acc","sens","spec","mcc","rmse","r")
regions <- c("East","Midwest","South")

## within region: 5-fold CV, trained only on the other folds of the same region
within <- matrix(NA,length(regions),7,dimnames=list(regions,MET))
for(R in regions){
  idxR <- which(reg==R & !is.na(y)); set.seed(7)
  fold <- sample(rep(1:5,length.out=length(idxR)))
  fm <- matrix(NA,5,7,dimnames=list(NULL,MET))
  for(k in 1:5){
    te <- idxR[fold==k]
    ym <- rep(NA_real_,n); tr <- idxR[fold!=k]; ym[tr] <- y[tr]   # same region only
    g <- gblup(ym); fm[k,] <- metvec(g[te], y[te])
  }
  within[R,] <- colMeans(fm,na.rm=TRUE)
  cat(sprintf("[within ] %-8s n=%3d auc=%.3f mcc=%+.3f sens=%.3f spec=%.3f r=%+.3f\n",
      R,length(idxR),within[R,"auc"],within[R,"mcc"],within[R,"sens"],within[R,"spec"],within[R,"r"]))
}

## across region: leave one region out (train on the other 2), without and with the PCs
across <- matrix(NA,length(regions),7,dimnames=list(regions,MET))
acrossPC <- matrix(NA,length(regions),7,dimnames=list(regions,MET))
for(R in regions){
  te <- which(reg==R & !is.na(y))
  ym <- y; ym[reg==R] <- NA
  across[R,]   <- metvec(gblup(ym,usePC=FALSE)[te], y[te])
  acrossPC[R,] <- metvec(gblup(ym,usePC=TRUE )[te], y[te])
  cat(sprintf("[across ] %-8s n=%3d auc=%.3f mcc=%+.3f sens=%.3f spec=%.3f r=%+.3f | +PC auc=%.3f mcc=%+.3f\n",
      R,length(te),across[R,"auc"],across[R,"mcc"],across[R,"sens"],across[R,"spec"],across[R,"r"],
      acrossPC[R,"auc"],acrossPC[R,"mcc"]))
}

## summary, mean (SD) over the 3 regions, and per region
fmt <- function(mat) sprintf("%.3f (%.3f)", colMeans(mat,na.rm=T), apply(mat,2,sd,na.rm=T))
summ <- data.frame(metric=MET,
                   within_region        = fmt(within),
                   across_region        = fmt(across),
                   across_region_plusPC = fmt(acrossPC))
cat("\n========== WITHIN vs ACROSS region (mean over 3 regions (SD)) ==========\n"); print(summ,row.names=FALSE)

perreg <- data.frame(region=regions,
  auc_within=round(within[,"auc"],3), auc_across=round(across[,"auc"],3), auc_acrossPC=round(acrossPC[,"auc"],3),
  mcc_within=round(within[,"mcc"],3), mcc_across=round(across[,"mcc"],3),
  sens_across=round(across[,"sens"],3), spec_across=round(across[,"spec"],3), r_across=round(across[,"r"],3))
cat("\n---- per region ----\n"); print(perreg,row.names=FALSE)

## write tables
write.csv(summ,  "reanalysis2022/data/within_vs_across_region_summary.csv", row.names=FALSE)
write.csv(perreg,"reanalysis2022/data/within_vs_across_region_perregion.csv", row.names=FALSE)
cat("\nsaved data/within_vs_across_region_summary.csv, within_vs_across_region_perregion.csv\n")
