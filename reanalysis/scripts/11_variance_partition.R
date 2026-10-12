# Step 3d: split pseudo-h2 into PRNP, region and polygenic remainder (rrBLUP mixed.solve, REML)
# PRNP G96S and region are fixed effects; the polygenic part uses the GRM without chr13 so PRNP is not counted twice

library(rrBLUP)
library(dplyr)

## working dir
ROOT <- "<dir>/Seabury_data/reanalysis"
S<-file.path(ROOT,"01_structure")

## read in GRMs (all autosomes and without chr13), eigenvalues floored at 1e-6
ids<-read.table(file.path(S,"grm.rel.id"))[,1]; n<-length(ids)
bendPSD<-function(K,eps=1e-6){ e<-eigen(K,symmetric=TRUE); v<-e$values; v[v<eps]<-eps
  K2<-e$vectors%*%diag(v)%*%t(e$vectors); (K2+t(K2))/2 }
Kall <-bendPSD(as.matrix(read.table(gzfile(file.path(S,"grm.rel.gz")))))
Kno13<-bendPSD(as.matrix(read.table(gzfile(file.path(S,"grm_no13.rel.gz")))))
meta<-read.delim(file.path(ROOT,"00_data/samples.tsv"),colClasses=c(clean_id="character")) %>% slice(match(ids,clean_id))
## PRNP genotypes (plink --recode A, see 00b_plink_2020.sh), missing set to the mean
G96<-read.table(file.path(ROOT, "work", "prnp.raw"),header=TRUE) %>% slice(match(ids,IID)) %>% pull(G96S_G_to_A_at286_A)
G96[is.na(G96)]<-mean(G96,na.rm=TRUE)
South<-as.numeric(meta$region=="South"); Mid<-as.numeric(meta$region=="Midwest")

## h2 under each model and variance partition of the full model
h2<-function(m) m$Vu/(m$Vu+m$Ve)
partition<-function(y,label){
  cat(sprintf("\n=================  %s  =================\n",label))
  base <-mixed.solve(y,K=Kall,  X=matrix(1,n,1))
  poly <-mixed.solve(y,K=Kno13, X=matrix(1,n,1))
  mp   <-mixed.solve(y,K=Kno13, X=cbind(1,G96))
  mpr  <-mixed.solve(y,K=Kno13, X=cbind(1,G96,South,Mid))
  cat(sprintf("h2 (GRM all autosomes, no fixed)         = %.3f   [reproduces paper ~0.64/0.59]\n",h2(base)))
  cat(sprintf("h2 (GRM ex-chr13, no fixed)              = %.3f\n",h2(poly)))
  cat(sprintf("h2 (GRM ex-chr13 | PRNP fixed)           = %.3f\n",h2(mp)))
  cat(sprintf("h2 (GRM ex-chr13 | PRNP+region fixed)    = %.3f\n",h2(mpr)))
  ## variance decomposition from the full model (PRNP + region fixed, polygenic GRM)
  b<-mpr$beta                       # intercept, PRNP, South, Mid
  fp<-G96*b[2]; fr<-South*b[3]+Mid*b[4]
  Vprnp<-var(fp); Vreg<-var(fr); Vcov<-2*cov(fp,fr); Vu<-mpr$Vu; Ve<-mpr$Ve
  Vp<-Vprnp+Vreg+Vcov+Vu+Ve
  cat(sprintf("\n  Variance decomposition (%% of phenotypic variance):\n"))
  cat(sprintf("    PRNP G96S (fixed)                 : %5.1f%%\n",100*Vprnp/Vp))
  cat(sprintf("    Region (fixed)                    : %5.1f%%\n",100*Vreg/Vp))
  cat(sprintf("    PRNP<->Region covariance          : %5.1f%%  (shared; they are correlated)\n",100*Vcov/Vp))
  cat(sprintf("    Polygenic remainder (GRM ex-chr13): %5.1f%%  <- the confounded 'heritability'\n",100*Vu/Vp))
  cat(sprintf("    Residual                          : %5.1f%%\n",100*Ve/Vp))
  cat(sprintf("    [check] var(y)=%.4f  model Vp=%.4f\n",var(y),Vp))
  data.frame(trait=label,PRNP=100*Vprnp/Vp,Region=100*Vreg/Vp,Cov=100*Vcov/Vp,Poly=100*Vu/Vp,Resid=100*Ve/Vp,
             h2_all=h2(base),h2_polyPRNPreg=h2(mpr))
}

## run for both traits
d1<-partition(meta$binary,"BINARY  CWD_Binary")
d2<-partition(meta$score, "SCORE   CMS_CWDScore")

## save
write.csv(bind_rows(d1,d2),file.path(S,"variance_partition.csv"),row.names=FALSE)
cat("\nsaved variance_partition.csv\n")
