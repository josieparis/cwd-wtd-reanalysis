# random 5-fold GBLUP in 2022 (10 replicates), then AUC with test animals binned by their maximum GRM
# relatedness to the training set (Fig 2b)
# writes reanalysis2022/data/relatedness_binned_auc.csv

library(rrBLUP)
library(dplyr)

## working dir
setwd("<dir>/Seabury_data")

source("reanalysis2022/scripts/load_grm.R")

## read in covariates and PLINK GRM
cov <- read.csv("reanalysis2022/data/cov_1170.csv", stringsAsFactors=FALSE)
K <- load_grm2022(cov$CMSID, bend = FALSE)                            # bent below after subsetting
y <- suppressWarnings(as.numeric(cov$binary))
keep <- !is.na(y); K <- K[keep,keep]; y <- y[keep]; n <- length(y)
e <- eigen(K, symmetric=TRUE); e$values[e$values < 1e-6] <- 1e-6      # bend to PSD
K <- e$vectors %*% diag(e$values) %*% t(e$vectors)
auc <- function(yy,pp){ r<-rank(pp); n1<-sum(yy==1); n0<-sum(yy==0); (sum(r[yy==1]) - n1*(n1+1)/2)/(n1*n0) }
Z <- diag(n)

## random 5-fold CV, 10 replicates (seed 1)
set.seed(1); ITER <- 10; KF <- 5
pred<-c(); yt<-c(); mr<-c()
for(it in 1:ITER){
  fold <- sample(rep(1:KF, length.out=n))
  for(k in 1:KF){
    te <- which(fold==k); tr <- which(fold!=k)
    sol <- mixed.solve(y=y[tr], Z=Z[tr,,drop=FALSE], K=K)
    pred <- c(pred, sol$u[te] + as.numeric(sol$beta)); yt <- c(yt, y[te])
    mr <- c(mr, apply(K[te,tr,drop=FALSE], 1, max))                   # max relatedness to the training set
  }
  cat(sprintf("iter %d done\n", it), file=stderr())
}

## AUC per relatedness bin
brk <- c(-Inf,0.05,0.125,0.25,Inf)
lab <- c("< 0.05","0.05-0.125","0.125-0.25","> 0.25")
deg <- c("unrelated","3rd-4th deg","2nd-3rd deg","1st-2nd deg")
out <- data.frame(bin=cut(mr, brk, labels=lab), yt, pred) %>%
  group_by(bin, .drop=FALSE) %>%                                       # keep empty bins
  summarise(n=n(), prev=round(mean(yt),3), auc=if(n>20) round(auc(yt,pred),3) else NA, .groups="drop") %>%
  mutate(bin=as.character(bin), degree=deg[match(bin, lab)], .after=bin) %>% as.data.frame()

## write table
write.csv(out, "reanalysis2022/data/relatedness_binned_auc.csv", row.names=FALSE)
print(out)
