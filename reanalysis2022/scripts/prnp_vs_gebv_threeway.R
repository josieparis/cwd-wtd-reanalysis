# leave-one-farm-out prediction for three predictors (Fig 3): PRNP codon 96 alone, the full GBLUP GEBV,
# and the GEBV from the GRM without chr13 (PRNP), i.e. the polygenic remainder
# writes reanalysis2022/data/prnp_threeway_summary.csv and prnp_threeway_perfarm.csv

library(rrBLUP)
library(dplyr)

## working dir
setwd("<dir>/Seabury_data")

source("reanalysis2022/scripts/load_grm.R")

## read in covariates and candidate genotypes
cov <- read.csv("reanalysis2022/data/cov_1170.csv", stringsAsFactors=FALSE) %>% mutate(y=suppressWarnings(as.numeric(binary)))
gt  <- read.delim("reanalysis2022/data/candidate_genos_1170.tsv", check.names=FALSE, stringsAsFactors=FALSE)
n <- nrow(cov); stopifnot(all(gt$CMSID==cov$CMSID))
y<-cov$y; farm<-cov$farm; ids<-cov$CMSID
prnp<-suppressWarnings(as.numeric(gt[["G96S_G_to_A_at286"]]))            # 0/1/2 protective allele dosage

## PLINK GRMs: all autosomes, and without chr13 (PRNP at 47,056,387)
Kfull <- load_grm2022(ids, "all"); Kno13 <- load_grm2022(ids, "no13")

## GEBV with kin.blup (mean + genomic value), and prediction metrics (cases called at >= 0.5)
gblup <- function(ymask,K){ fit<-kin.blup(data.frame(gid=ids,yy=ymask),geno="gid",pheno="yy",K=K)
                            g<-fit$g[match(ids,names(fit$g))]; mean(ymask,na.rm=TRUE)+g }
metvec <- function(pred,truth,thr=0.5){ ok<-!is.na(pred)&!is.na(truth); p<-pred[ok]; t<-truth[ok]
  np<-sum(t==1); nn<-sum(t==0); N<-length(t); if(!np||!nn) return(c(auc=NA,acc=NA,sens=NA,spec=NA,mcc=NA,rmse=NA,r=NA))
  auc<-(sum(rank(p)[t==1])-np*(np+1)/2)/(np*nn); rmse<-sqrt(mean((p-t)^2)); r<-if(sd(p)>0)cor(p,t) else NA_real_
  cls<-as.integer(p>=thr); TP<-sum(cls==1&t==1);TN<-sum(cls==0&t==0);FP<-sum(cls==1&t==0);FN<-sum(cls==0&t==1)
  acc<-(TP+TN)/N; sens<-TP/(TP+FN); spec<-TN/(TN+FP); den<-sqrt(as.double(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  c(auc=auc,acc=acc,sens=sens,spec=spec,mcc=if(den>0)(TP*TN-FP*FN)/den else 0,rmse=rmse,r=r) }
MET<-c("auc","acc","sens","spec","mcc","rmse","r")

## leave one farm out (farms with >= 10 deer and >= 3 cases and controls), with deer and prevalence per farm
info<-cov %>% filter(!is.na(y), !is.na(farm), farm!="") %>% group_by(farm) %>%
  summarise(n=n(), cases=sum(y==1), controls=sum(y==0), prev=mean(y==1), .groups="drop") %>%
  filter(n>=10, cases>=3, controls>=3) %>% select(farm, n, prev) %>% as.data.frame()
qfarms<-info$farm; cat("[chk] qualifying farms:",length(qfarms),"\n\n")
mods<-c("PRNP_only","GEBV_full","GEBV_minusPRNPchr")
pf<-lapply(mods,function(m) matrix(NA,length(qfarms),7,dimnames=list(qfarms,MET))); names(pf)<-mods
for(fi in seq_along(qfarms)){ F<-qfarms[fi]; fidx<-which(farm==F&!is.na(y))
  ym<-y; ym[farm==F]<-NA
  pf$GEBV_full[fi,]         <- metvec(gblup(ym,Kfull)[fidx], y[fidx])
  pf$GEBV_minusPRNPchr[fi,] <- metvec(gblup(ym,Kno13)[fidx], y[fidx])
  # PRNP alone: logistic regression fitted on the other farms
  tr<-which(farm!=F&!is.na(y)&!is.na(prnp)); fit<-glm(yy~d,binomial,data.frame(yy=y[tr],d=prnp[tr]))
  pf$PRNP_only[fi,]         <- metvec(as.numeric(predict(fit,data.frame(d=prnp[fidx]),type="response")), y[fidx])
  cat(sprintf("  %-9s n=%3d prev=%.2f | PRNP %.3f  GEBV_full %.3f  GEBV_minusPRNPchr %.3f\n",
      F,info$n[fi],info$prev[fi],pf$PRNP_only[fi,"auc"],pf$GEBV_full[fi,"auc"],pf$GEBV_minusPRNPchr[fi,"auc"])) }

## summary over farms, mean (SD), and AUC per farm
fmt<-function(m) sprintf("%.3f (%.3f)",colMeans(m,na.rm=T),apply(m,2,sd,na.rm=T))
summ<-data.frame(metric=MET, PRNP_only=fmt(pf$PRNP_only), GEBV_full=fmt(pf$GEBV_full), GEBV_minusPRNPchr=fmt(pf$GEBV_minusPRNPchr))
cat("\n=========== THREE-WAY out-of-farm (mean over 13 farms (SD)) ===========\n"); print(summ,row.names=FALSE)
out<-info %>% mutate(auc_PRNP=round(pf$PRNP_only[,"auc"],3), auc_GEBVfull=round(pf$GEBV_full[,"auc"],3),
                     auc_GEBVminus13=round(pf$GEBV_minusPRNPchr[,"auc"],3))
cat("\n---- per farm (out-of-farm AUC) ----\n"); print(out,row.names=FALSE)

## write tables
write.csv(summ,"reanalysis2022/data/prnp_threeway_summary.csv",row.names=FALSE)
write.csv(out,"reanalysis2022/data/prnp_threeway_perfarm.csv",row.names=FALSE)
cat("\nsaved reanalysis2022/data/prnp_threeway_summary.csv, prnp_threeway_perfarm.csv\n")
