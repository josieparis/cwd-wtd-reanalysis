# GxEMM (Dahl et al. 2020) with region as the environment, with and without farm as a fixed effect (2022 CWD data)
# Model A: region + sex + age (as in run_gxemm_cluster.R), Model B: farm + sex + age, on the same deer
# robustness check: 99.7% of first-degree pairs share a farm, so the farm-adjusted h2 is an upper bound

library(GxEMM)

## run inside the cluster_gxemm folder (24 GxEMM fits, ask for 6-8 h)
LDAK <- "./ldak5.linux "                       # trailing space needed (GxEMM appends args)
stopifnot(file.exists("ldak5.linux"))

## read in data
n   <- 1171
K   <- matrix(readBin("K_1170.bin","double",n*n), n, n, byrow=TRUE)   # PLINK GRM, from build_K_from_plink.R
cov <- read.csv("cov_1170.csv", stringsAsFactors=FALSE)
cov$cwd  <- suppressWarnings(as.numeric(cov$binary))
cov$sexn <- suppressWarnings(as.numeric(cov$sex))

## complete cases on farms with >= 10 deer and both cases and controls
tb   <- table(cov$farm[!is.na(cov$cwd)])
both <- tapply(cov$cwd[!is.na(cov$cwd)], cov$farm[!is.na(cov$cwd)], function(v) length(unique(v))>1)
keepfarm <- names(tb)[ tb>=10 & both[names(tb)] ]
ok <- !is.na(cov$cwd) & !is.na(cov$sexn) & !is.na(cov$age) & cov$farm %in% keepfarm
K  <- K[ok,ok]; cov <- cov[ok,]; N <- nrow(K)

y    <- cov$cwd
reg  <- factor(cov$region, levels=c("East","Midwest","South")); lv <- levels(reg)
farm <- factor(cov$farm)
Z    <- model.matrix(~0+reg); colnames(Z) <- lv                 # N x 3 environment indicators

XA <- cbind( Z[,-1], sex=cov$sexn, age=cov$age )                # Model A: region-adjusted (reference)
XB <- model.matrix(~ farm + sexn + age, data=cov)[,-1,drop=FALSE] # Model B: farm-adjusted (subsumes region)

cat(sprintf("N=%d across %d farms, %d regions\n", N, nlevels(farm), ncol(Z)))
cat(sprintf("farms kept: %s\n", paste(keepfarm, collapse=", ")))
cat(sprintf("overall prevalence=%.3f | by region: %s\n", mean(y),
    paste(sprintf("%s=%.2f",lv,tapply(y,reg,mean)),collapse=" ")))
cat(sprintf("design ranks: A=%d/%d cols, B=%d/%d cols\n\n",
    qr(cbind(1,XA))$rank, ncol(XA)+1, qr(cbind(1,XB))$rank, ncol(XB)+1))
stopifnot(qr(cbind(1,XB))$rank == ncol(XB)+1)   # farm design must be full rank

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
  invisible(p)
}

## fit HOM, IID and FREE on the observed scale and at three liability prevalences
runset <- function(tag, X){
  res <- list()
  oh <- GxEMM(y, X, K, Z, gtype='hom',                ldak_loc=LDAK)
  oi <- GxEMM(y, X, K, Z, gtype='iid',                ldak_loc=LDAK)
  of <- GxEMM(y, X, K, Z, gtype='free', etype='free', ldak_loc=LDAK)
  report(sprintf("%s | OBSERVED SCALE", tag), oh, oi, of)
  res$obs <- list(hom=oh, iid=oi, free=of)

  ycc <- y + 1L    # binary=TRUE wants 1 = control, 2 = case
  for(prev in c(0.05, 0.10, 0.26)){
    lh <- GxEMM(ycc, X, K, Z, gtype='hom',                binary=TRUE, prev=prev, ldak_loc=LDAK)
    li <- GxEMM(ycc, X, K, Z, gtype='iid',                binary=TRUE, prev=prev, ldak_loc=LDAK)
    lf <- GxEMM(ycc, X, K, Z, gtype='free', etype='free', binary=TRUE, prev=prev, ldak_loc=LDAK)
    report(sprintf("%s | LIABILITY SCALE (prev=%.2f)", tag, prev), lh, li, lf)
    res[[sprintf("liab_%.2f",prev)]] <- list(hom=lh, iid=li, free=lf)
  }
  res
}

## run models and save
out <- list()
out$A_region_adj <- runset("MODEL A (region-adjusted, reference)", XA)
out$B_farm_adj   <- runset("MODEL B (FARM-ADJUSTED)",              XB)
out$meta <- list(N=N, farms=keepfarm, prev=mean(y), prev_by_region=tapply(y,reg,mean))
saveRDS(out, "gxemm_farmadj_results.rds")

cat("\n\n=================== WHAT TO LOOK FOR ===================\n")
cat("1) HOM liability h2 under MODEL B should fall into a PLAUSIBLE range (<1) where Model A gave >1.\n")
cat("   If so, the impossible values were herd/transmission variance being forced into the genetic term.\n")
cat("2) The FREE-model GxE p-value on the LIABILITY scale should REMAIN NON-SIGNIFICANT under Model B.\n")
cat("   That is the robustness result: no gene x region interaction, with or without farm adjustment.\n")
cat("3) If Model B's h2 collapses toward 0 with huge SEs, that is the near-collinearity of farm and\n")
cat("   kinship (99.7%% of 1st-degree pairs are co-resident) -- report as NON-IDENTIFIABLE, not as 0.\n")
cat("\nSaved gxemm_farmadj_results.rds\n")
