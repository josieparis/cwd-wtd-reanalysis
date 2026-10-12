# GxEMM (Dahl et al. 2020) on the 2022 CWD data with region as the environment,
# on the observed and the liability scale (binary = TRUE), using the bundled ldak5.linux
# first-time setup in README_CLUSTER.md (install msm and GxEMM, make ldak5.linux executable)

library(GxEMM)

## run inside the cluster_gxemm folder
LDAK <- "./ldak5.linux "                       # trailing space needed (GxEMM appends args)
stopifnot(file.exists("ldak5.linux"))

## read in data
n   <- 1171
K   <- matrix(readBin("K_1170.bin","double",n*n), n, n, byrow=TRUE)   # PLINK GRM, from build_K_from_plink.R
cov <- read.csv("cov_1170.csv")
cov$cwd  <- suppressWarnings(as.numeric(cov$binary))
cov$sexn <- suppressWarnings(as.numeric(cov$sex))
ok  <- !is.na(cov$cwd) & !is.na(cov$sexn) & !is.na(cov$age)
K   <- K[ok,ok]; cov <- cov[ok,]; N <- nrow(K)
y   <- cov$cwd
reg <- factor(cov$region, levels=c("East","Midwest","South")); lv <- levels(reg)
Z   <- model.matrix(~0+reg); colnames(Z) <- lv          # N x 3 environment indicators
X   <- cbind( Z[,-1], sex=cov$sexn, age=cov$age )        # Z (minus 1 column) + covariates
cat(sprintf("N=%d  regions: %s\n", N, paste(colnames(Z),colSums(Z),collapse=" ")))
cat(sprintf("sample prevalence overall=%.3f  by region: %s\n\n", mean(y),
    paste(sprintf("%s=%.2f",lv,tapply(y,reg,mean)),collapse=" ")))

## print the HOM, IID and FREE results
report <- function(tag, out_hom, out_iid, out_free){
  cat("\n==================", tag, "==================\n")
  cat(sprintf("HOM  h2 = %.3f (SE %.3f)  Wald p(h2>0) = %.2e\n",
      out_hom$h2[1], sqrt(out_hom$h2Covmat[1,1]), Waldtest(out_hom$h2[1], out_hom$h2Covmat[1,1])))
  cat(sprintf("IID  shared-GxE h2 = %.3f (SE %.3f)   Wald p(GxE) = %.3g\n",
      out_iid$h2[2], sqrt(out_iid$h2Covmat[2,2]), Waldtest(out_iid$h2[2], out_iid$h2Covmat[2,2])))
  enms <- gsub("East","Northeast",colnames(Z))
  fse  <- sqrt(pmax(diag(as.matrix(out_free$h2Covmat)), 0))
  cat("FREE per-environment h2 (estimate +/- SE):\n")
  for(i in seq_along(out_free$h2))
    cat(sprintf("     %-12s h2 = %8.3f  SE = %8.3f\n", enms[i], out_free$h2[i], fse[i]))
  ng <- ncol(Z); het <- 2:(1+ng)
  cat("FREE sig2s:\n"); print(round(out_free$sig2s,4))
  p <- tryCatch(MVWaldtest(out_free$sig2s[het], out_free$sig2Var[het,het]), error=function(e) NA)
  cat(sprintf("FREE  MV-Wald p for genetic heterogeneity (GxE) = %.3g\n", p))
}

## observed scale
oh <- GxEMM(y, X, K, Z, gtype='hom',                ldak_loc=LDAK)
oi <- GxEMM(y, X, K, Z, gtype='iid',                ldak_loc=LDAK)
of <- GxEMM(y, X, K, Z, gtype='free', etype='free', ldak_loc=LDAK)
report("OBSERVED SCALE", oh, oi, of)
res <- list(obs=list(hom=oh,iid=oi,free=of))

## liability scale (PCGC) at three population prevalences
ycc <- y + 1L    # binary=TRUE wants 1 = control, 2 = case
for(prev in c(0.05, 0.10, 0.26)){          # 0.26 = their weighted mean farm prevalence
  lh <- GxEMM(ycc, X, K, Z, gtype='hom',                binary=TRUE, prev=prev, ldak_loc=LDAK)
  li <- GxEMM(ycc, X, K, Z, gtype='iid',                binary=TRUE, prev=prev, ldak_loc=LDAK)
  lf <- GxEMM(ycc, X, K, Z, gtype='free', etype='free', binary=TRUE, prev=prev, ldak_loc=LDAK)
  report(sprintf("LIABILITY SCALE (prev=%.2f)",prev), lh, li, lf)
  res[[sprintf("liab_%.2f",prev)]] <- list(hom=lh,iid=li,free=lf)
}

## save
saveRDS(res, "gxemm_cluster_results.rds")
cat("\nSaved gxemm_cluster_results.rds\n")
cat("\nINTERPRETATION: if the FREE-model GxE p-value is NON-significant on the LIABILITY scale,\n")
cat("the region 'GxE' seen on the observed scale was the binary-prevalence artefact (+ farm structure),\n")
cat("i.e. no real gene x region interaction — consistent with the local analysis.\n")
