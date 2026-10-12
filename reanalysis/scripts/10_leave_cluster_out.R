# Step 3c: leave-one-genetic-cluster-out GBLUP (in place of leave-one-herd-out)
# clusters = connected components of the GRM cut at t; whole clusters go to one of 5 folds, so no relative
# at >= t of a test animal is in training. Lower t removes more distant relatives; compare AUC with random CV

library(rrBLUP)
library(pROC)
library(dplyr)

## working dir
ROOT <- "<dir>/Seabury_data/reanalysis"
S<-file.path(ROOT,"01_structure")
auc<-function(y,p) as.numeric(pROC::auc(response=y,predictor=p,quiet=TRUE,levels=c(0,1),direction="<"))

## read in GRM and phenotypes
ids<-read.table(file.path(S,"grm.rel.id"))[,1]; n<-length(ids)
K<-as.matrix(read.table(gzfile(file.path(S,"grm.rel.gz")))); dimnames(K)<-list(ids,ids)
Kd<-K+diag(1e-4,n); Koff<-K; diag(Koff)<- -Inf
meta<-read.delim(file.path(ROOT,"00_data/samples.tsv"),colClasses=c(clean_id="character")) %>% slice(match(ids,clean_id))
y<-meta$binary

## GBLUP prediction of the test animals
gblup_pred<-function(tr,te){ df<-data.frame(gid=ids,yv=y) %>% mutate(yv=replace(yv,te,NA))
  kin.blup(data=df,geno="gid",pheno="yv",K=Kd)$g[te] }

## clusters = connected components of GRM >= thr
components<-function(thr){
  parent<-1:n; find<-function(x){ while(parent[x]!=x) x<-parent[x]; x }
  idx<-which(K>=thr & upper.tri(K),arr.ind=TRUE)
  for(r in seq_len(nrow(idx))){ a<-idx[r,1];b<-idx[r,2];ra<-find(a);rb<-find(b); if(ra!=rb) parent[rb]<-ra }
  as.integer(factor(sapply(1:n,find)))
}
## whole clusters to F folds, in random order, each to the fold with fewest animals
assign_folds<-function(clust,F,seed){ set.seed(seed); nc<-max(clust); sz<-tabulate(clust,nc)
  ord<-sample(nc); load<-rep(0,F); cf<-integer(nc)
  for(c in ord){ f<-which.min(load); cf[c]<-f; load[f]<-load[f]+sz[c] }
  cf[clust] }

## CV over the given folds: AUC and mean max relatedness of test animals to training
blocked_cv<-function(folds){ pred<-rep(NA,n); mr<-rep(NA,n)
  for(k in sort(unique(folds))){ te<-which(folds==k); tr<-which(folds!=k)
    pred[te]<-gblup_pred(tr,te); mr[te]<-apply(Koff[te,tr,drop=FALSE],1,max) }
  c(auc=auc(y,pred), mean_maxrel=mean(mr)) }

## run (8 seeds per scheme)
SEEDS<-1:8
cat("=== Leave-one-genetic-cluster-out GBLUP (binary) ===\n")
cat(sprintf("%-22s %8s %8s %10s %12s\n","scheme","AUC","AUC_sd","nclusters","maxrel2train"))
## random 5-fold baseline
rb<-t(sapply(SEEDS,function(s){set.seed(100+s); f<-sample(rep(1:5,length.out=n)); blocked_cv(f)}))
cat(sprintf("%-22s %8.3f %8.3f %10s %12.3f\n","random 5-fold",mean(rb[,1]),sd(rb[,1]),"-",mean(rb[,2])))
## cluster-out CV at t = 0.45, 0.35 and 0.25
for(thr in c(0.45,0.35,0.25)){
  cl<-components(thr); nc<-max(cl)
  res<-t(sapply(SEEDS,function(s){ f<-assign_folds(cl,5,s); blocked_cv(f) }))
  cat(sprintf("cluster-out t>=%-4.2f      %8.3f %8.3f %10d %12.3f\n",thr,mean(res[,1]),sd(res[,1]),nc,mean(res[,2])))
}
cat("\nRef: Test B fully-unrelated stratum (GRM<0.125 to training) AUC=0.48 (chance)\n")
cat("Interpretation: as t drops, more relatives removed from training (maxrel2train falls) and AUC falls toward chance.\n")
