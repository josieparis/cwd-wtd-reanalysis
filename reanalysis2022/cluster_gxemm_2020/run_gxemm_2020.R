# GxEMM PCGC on the 2020 dataset (n = 807), region as the environment and region + sex + age as covariates
# 2020 has no farm IDs, so only the region-adjusted model; the HOM liability h2 is the 2020 PCGC row of Table 4

library(GxEMM)

## run inside the cluster_gxemm_2020 folder
LDAK <- "./ldak5.linux "
stopifnot(file.exists("ldak5.linux"))

## read in data
n   <- 807
K   <- matrix(readBin("K_807.bin","double",n*n), n, n, byrow=TRUE)   # PLINK GRM, from build_K_from_plink.R (positional align)
cov <- read.csv("cov_807.csv")
cov$cwd  <- suppressWarnings(as.numeric(cov$binary))
cov$sexn <- suppressWarnings(as.numeric(cov$sex))
ok  <- !is.na(cov$cwd) & !is.na(cov$sexn) & !is.na(cov$age)
K   <- K[ok,ok]; cov <- cov[ok,]; N <- nrow(K)
y   <- cov$cwd
reg <- factor(cov$region, levels=c("East","Midwest","South")); lv <- levels(reg)
Z   <- model.matrix(~0+reg); colnames(Z) <- lv
X   <- cbind( Z[,-1], sex=cov$sexn, age=cov$age )
cat(sprintf("N=%d  regions: %s\n", N, paste(colnames(Z),colSums(Z),collapse=" ")))
cat(sprintf("sample prevalence overall=%.3f  by region: %s\n\n", mean(y),
    paste(sprintf("%s=%.2f",lv,tapply(y,reg,mean)),collapse=" ")))

## HOM h2 and SE
hom_h2 <- function(o) c(h2=o$h2[1], se=sqrt(o$h2Covmat[1,1]))
res <- list()

## observed-scale HOM (reference)
oh <- GxEMM(y, X, K, Z, gtype='hom', ldak_loc=LDAK)
res$obs_hom <- oh
cat(sprintf("OBSERVED HOM h2 = %.3f (SE %.3f)\n\n", hom_h2(oh)[1], hom_h2(oh)[2]))

## liability-scale PCGC (region-adjusted) at each prevalence, the Table 4 row
ycc <- y + 1L
cat("=== region-adjusted PCGC (HOM, liability) ===\n")
for(prev in c(0.05, 0.10, 0.26)){
  lh <- GxEMM(ycc, X, K, Z, gtype='hom', binary=TRUE, prev=prev, ldak_loc=LDAK)
  res[[sprintf("liab_hom_%.2f",prev)]] <- lh
  cat(sprintf("prev=%.2f : h2_L = %.3f (SE %.3f)\n", prev, hom_h2(lh)[1], hom_h2(lh)[2]))
}

## region GxE test (FREE model, liability, prev 0.05), can be unstable (East n = 35) so wrapped in tryCatch
cat("\n=== region GxE test (liability FREE, prev=0.05) ===\n")
tryCatch({
  lf <- GxEMM(ycc, X, K, Z, gtype='free', etype='free', binary=TRUE, prev=0.05, ldak_loc=LDAK)
  res$liab_free_0.05 <- lf
  ng <- ncol(Z); het <- 2:(1+ng)
  cat("FREE region-specific h2 (estimate +/- SE):\n")
  enms <- gsub("East","Northeast",colnames(Z))
  fse  <- sqrt(pmax(diag(as.matrix(lf$h2Covmat)), 0))
  for(i in seq_along(lf$h2))
    cat(sprintf("     %-12s h2 = %8.3f  SE = %8.3f\n", enms[i], lf$h2[i], fse[i]))
  p <- tryCatch(MVWaldtest(lf$sig2s[het], lf$sig2Var[het,het]), error=function(e) NA)
  cat(sprintf("FREE MV-Wald p for gene x region (GxE) = %.3g\n", p))
}, error=function(e) cat("FREE model failed:", conditionMessage(e), "\n"))

## save
saveRDS(res, "gxemm_2020_results.rds")
cat("\nSaved gxemm_2020_results.rds\n")
