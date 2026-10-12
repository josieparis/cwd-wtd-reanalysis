# Reproduce Seabury Table 1: GBLUP prediction with k = 3 and 5 folds, no covariates vs sex + age + region,
# 50 iterations, classification threshold 0.5 (rrBLUP mixed.solve)
# test predictions are X_test %*% beta + g_test, i.e. fixed effects added back; metrics are mean (SD)

library(rrBLUP)
library(pROC)

## working dir
ROOT <- "<dir>/Seabury_data/reanalysis"
S<-file.path(ROOT,"01_structure")

## read in GRM and phenotypes
ids<-read.table(file.path(S,"grm.rel.id"))[,1]; n<-length(ids)
K<-as.matrix(read.table(gzfile(file.path(S,"grm.rel.gz")))); dimnames(K)<-list(ids,ids); K<-K+diag(1e-4,n)
meta<-read.delim(file.path(ROOT,"00_data/samples.tsv"),colClasses=c(clean_id="character")); rownames(meta)<-meta$clean_id; meta<-meta[ids,]
y<-meta$binary

## fixed effects: intercept only, or intercept + sex + age + region
Xnone<-matrix(1,n,1)
Xcov <-model.matrix(~ factor(sex_raw)+as.numeric(age)+factor(region), data=meta)
In<-diag(n)
THR<-0.5; ITER<-50

## predict a test fold: fixed effects + breeding value
pred_fold<-function(train,test,X){
  Ztr<-In[train,,drop=FALSE]                                  # incidence, n_train x n
  fit<-mixed.solve(y=y[train], Z=Ztr, K=K, X=X[train,,drop=FALSE])
  as.numeric(X[test,,drop=FALSE]%*%as.numeric(fit$beta)) + as.numeric(fit$u)[test]
}
## one CV iteration: AUC, MCC, accuracy, sensitivity, specificity, RMSE
one_iter<-function(Kfold,X){
  fold<-sample(rep(1:Kfold,length.out=n)); pred<-rep(NA,n)
  for(k in 1:Kfold){ te<-which(fold==k); pred[te]<-pred_fold(which(fold!=k),te,X) }
  yh<-as.integer(pred>THR)
  TP<-sum(yh&y);TN<-sum(!yh&!y);FP<-sum(yh&!y);FN<-sum(!yh&y)
  sens<-TP/(TP+FN); spec<-TN/(TN+FP); acc<-(TP+TN)/n
  den<-sqrt(as.double(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)); mcc<-ifelse(den==0,0,(TP*TN-FP*FN)/den)
  c(AUC=as.numeric(pROC::auc(y,pred,quiet=TRUE,levels=c(0,1),direction="<")),
    MCC=mcc,Accuracy=acc,Sensitivity=sens,Specificity=spec,RMSE=sqrt(mean((pred-y)^2)))
}

## run the four configurations
set.seed(20200304)
cfgs<-list(list(3,Xnone,"k=3  None"),list(3,Xcov,"k=3  Sex,Age,Region"),
           list(5,Xnone,"k=5  None"),list(5,Xcov,"k=5  Sex,Age,Region"))
cat(sprintf("%-22s %14s %14s %14s %14s %14s %14s\n","config","AUC","Matthews","Accuracy","Sensitivity","Specificity","RMSE"))
out<-list()
for(cf in cfgs){
  M<-t(replicate(ITER, one_iter(cf[[1]],cf[[2]]))); m<-colMeans(M); s<-apply(M,2,sd)
  cat(sprintf("%-22s %6.4f(%.4f) %6.4f(%.4f) %6.4f(%.4f) %6.4f(%.4f) %6.4f(%.4f) %6.4f(%.4f)\n",
      cf[[3]],m[1],s[1],m[2],s[2],m[3],s[3],m[4],s[4],m[5],s[5],m[6],s[6]))
  out[[cf[[3]]]]<-data.frame(config=cf[[3]],metric=names(m),mean=m,sd=s)
}

## save
write.csv(do.call(rbind,out),file.path(S,"table1_reproduction.csv"),row.names=FALSE)

## values reported in Seabury Table 1
cat("\nSeabury Table 1 (reported):\n")
cat("k=3 None : AUC .8471  MCC .5871  Acc .8167  Sens .6447  Spec .9101  RMSE .3768\n")
cat("k=3 Cov  : AUC .8534  MCC .5787  Acc .8119  Sens .6735  Spec .8870  RMSE .3746\n")
cat("k=5 None : AUC .8485  MCC .5940  Acc .8198  Sens .6496  Spec .9121  RMSE .3754\n")
cat("k=5 Cov  : AUC .8542  MCC .5870  Acc .8159  Sens .6716  Spec .8942  RMSE .3734\n")
