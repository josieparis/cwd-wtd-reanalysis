# GxEMM (Dahl et al. 2020) with farm as the environment instead of region (2022 CWD data)
# farms with >= 10 deer and both cases and controls; the FREE model (~2 x Nfarm + 1 variance components)
# is unstable or does not converge, so per-farm h2 are not reliable and the IID GxE test is the one to use

library(GxEMM)

## run inside the cluster_gxemm folder
LDAK <- "./ldak5.linux "
stopifnot(file.exists("ldak5.linux"))

## read in data
n <- 1171
K <- matrix(readBin("K_1170.bin","double",n*n), n, n, byrow=TRUE)   # PLINK GRM, from build_K_from_plink.R
cov <- read.csv("cov_1170.csv")
cov$cwd  <- suppressWarnings(as.numeric(cov$binary)); cov$sexn <- suppressWarnings(as.numeric(cov$sex))
ok <- !is.na(cov$cwd) & !is.na(cov$sexn) & !is.na(cov$age)
K <- K[ok,ok]; cov <- cov[ok,]

## farms with >= 10 deer and both cases and controls (h2 is undefined at 0% or 100% CWD)
ft <- table(cov$farm); prev <- tapply(cov$cwd, cov$farm, mean)
keep_farm <- names(ft)[ ft>=10 & prev>0 & prev<1 ]
cat("Farms (>=10 deer):\n"); print(round(sort(prev[names(ft)[ft>=10]]),2))
cat("\nEXCLUDED (0% or 100% CWD, h2 undefined):",
    paste(setdiff(names(ft)[ft>=10], keep_farm), collapse=", "), "\n")
sel <- cov$farm %in% keep_farm
K <- K[sel,sel]; cov <- cov[sel,]; N <- nrow(K)
cat(sprintf("\nAnalysing %d deer across %d farms\n\n", N, length(keep_farm)))

y   <- cov$cwd
frm <- factor(cov$farm, levels=keep_farm)
Z   <- model.matrix(~0+frm); colnames(Z) <- keep_farm     # N x Nfarm environment indicators
X   <- cbind( Z[,-1], sex=cov$sexn, age=cov$age )          # Z (minus 1 column) + covariates

## GxEMM returns the string "LDAK did not converge" instead of a list when REML fails
## the FREE farm model is not expected to converge (99.7% of first-degree relatives are herd-mates), so keep going
ok_fit <- function(o) is.list(o) && !is.null(o$h2)
safe   <- function(expr) tryCatch(expr, error=function(e) paste("ERROR:", conditionMessage(e)))

## print the HOM, IID and FREE results
report <- function(tag, out_hom, out_iid, out_free){
  cat("\n==================", tag, "==================\n")
  if(ok_fit(out_hom)) cat(sprintf("HOM  h2 = %.3f (SE %.3f)  Wald p(h2>0) = %.2e\n",
      out_hom$h2[1], sqrt(out_hom$h2Covmat[1,1]), Waldtest(out_hom$h2[1], out_hom$h2Covmat[1,1])))
  else cat("HOM  -- DID NOT CONVERGE --\n")
  if(ok_fit(out_iid)) cat(sprintf("IID  shared-GxE h2 = %.3f (SE %.3f)   Wald p(GxE) = %.3g\n",
      out_iid$h2[2], sqrt(out_iid$h2Covmat[2,2]), Waldtest(out_iid$h2[2], out_iid$h2Covmat[2,2])))
  else cat("IID  -- DID NOT CONVERGE --\n")
  if(ok_fit(out_free)){
    ng <- ncol(Z); het <- 2:(1+ng)
    p <- tryCatch(MVWaldtest(out_free$sig2s[het], out_free$sig2Var[het,het]), error=function(e) NA)
    cat(sprintf("FREE per-farm h2 range: [%.2f, %.2f] (unstable; see caveats)\n",
        suppressWarnings(min(out_free$h2)), suppressWarnings(max(out_free$h2))))
    { enms <- gsub("East","Northeast",colnames(Z))
      fse  <- sqrt(pmax(diag(as.matrix(out_free$h2Covmat)), 0))
      cat("FREE per-farm h2 (estimate +/- SE):\n")
      for(i in seq_along(out_free$h2))
        cat(sprintf("     %-14s h2 = %8.3f  SE = %8.3f\n", enms[i], out_free$h2[i], fse[i])) }
    cat(sprintf("FREE  MV-Wald p for genetic heterogeneity ACROSS FARMS (GxE) = %.3g\n", p))
  } else {
    cat("FREE -- DID NOT CONVERGE (expected: ~2*Nfarm+1 variance components are not identifiable\n")
    cat("        when farms and families coincide). Report as NON-IDENTIFIABLE, and rely on the IID\n")
    cat("        shared-GxE test above, which is estimable.\n")
  }
}

res <- list(farms=keep_farm, N=N)

## observed scale
oh<-safe(GxEMM(y,X,K,Z,gtype='hom',ldak_loc=LDAK)); oi<-safe(GxEMM(y,X,K,Z,gtype='iid',ldak_loc=LDAK))
of<-safe(GxEMM(y,X,K,Z,gtype='free',etype='free',ldak_loc=LDAK))
report("OBSERVED SCALE (farm environment)", oh, oi, of)
res$obs <- list(hom=oh,iid=oi,free=of)
saveRDS(res, "gxemm_farm_results.rds")            # save as we go

## liability scale (PCGC)
ycc <- y + 1L
for(prev in c(0.05, 0.10, 0.26)){
  lh<-safe(GxEMM(ycc,X,K,Z,gtype='hom',binary=TRUE,prev=prev,ldak_loc=LDAK))
  li<-safe(GxEMM(ycc,X,K,Z,gtype='iid',binary=TRUE,prev=prev,ldak_loc=LDAK))
  lf<-safe(GxEMM(ycc,X,K,Z,gtype='free',etype='free',binary=TRUE,prev=prev,ldak_loc=LDAK))
  report(sprintf("LIABILITY SCALE prev=%.2f (farm environment)",prev), lh, li, lf)
  res[[sprintf("liab_%.2f",prev)]] <- list(hom=lh,iid=li,free=lf)
  saveRDS(res, "gxemm_farm_results.rds")
}
cat("\nSaved gxemm_farm_results.rds\n")
cat("INTERPRETATION: the estimable test here is the IID shared-GxE p-value. Non-significant => no\n")
cat("farm-level GxE (herd differences are allele-frequency/prevalence/relatedness structure, not gene x herd).\n")
cat("FREE non-convergence is itself a finding: per-farm genetic variances are NOT identifiable in this design.\n")
