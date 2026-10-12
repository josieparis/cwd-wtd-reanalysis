# heritability with standard errors (rrBLUP, PLINK GRM) for every Table 3 model in both cohorts
# rrBLUP gives no SEs, so the SE comes from the REML information matrix and the delta method
# writes reanalysis2022/cluster_h2/rrblup/h2_rrblup_se.csv

library(rrBLUP)
library(dplyr)

## working dir
setwd("<dir>/Seabury_data")

OUT <- "reanalysis2022/cluster_h2/rrblup"; dir.create(OUT, showWarnings = FALSE)

## h2 and its SE (expected information matrix, then delta method for h2 = s2g/(s2g+s2e))
h2_se <- function(y, K, X){
  ok <- !is.na(y) & stats::complete.cases(X)
  y <- y[ok]; K <- K[ok, ok, drop=FALSE]; X <- X[ok, , drop=FALSE]
  ms <- mixed.solve(y, K=K, X=X)
  s2g <- as.numeric(ms$Vu); s2e <- as.numeric(ms$Ve); n <- length(y)
  V    <- s2g*K + s2e*diag(n)
  Vi   <- chol2inv(chol(V))
  XtVi <- crossprod(X, Vi)
  P    <- Vi - t(XtVi) %*% solve(XtVi %*% X) %*% XtVi
  PK   <- P %*% K
  I11  <- 0.5*sum(PK * t(PK))          # tr(PK PK)
  I12  <- 0.5*sum(PK * t(P))           # tr(PK P)
  I22  <- 0.5*sum(P  * t(P))           # tr(P P)
  Vth  <- solve(matrix(c(I11,I12,I12,I22), 2, 2))
  tot  <- s2g + s2e
  J    <- c(s2e/tot^2, -s2g/tot^2)     # d h2 / d(s2g, s2e)
  c(h2 = s2g/tot, se = sqrt(as.numeric(t(J) %*% Vth %*% J)), n = n)
}
read_rel <- function(prefix){                       # PLINK square GRM, dimnames = IID
  K <- as.matrix(read.table(gzfile(paste0(prefix, ".rel.gz"))))
  ids <- read.table(paste0(prefix, ".rel.id"), stringsAsFactors = FALSE)$V2
  dimnames(K) <- list(ids, ids); K
}
pd <- function(K){                                  # floor eigenvalues so K is positive definite
  ei <- eigen(K, symmetric=TRUE); ev <- pmax(ei$values, 1e-6)
  ei$vectors %*% (ev * t(ei$vectors))
}
## print and store each result
res <- list()
add <- function(cohort, trait, covset, r){
  cat(sprintf("  %s %-7s %-38s h2 = %.3f +- %.3f   (n=%d)\n", cohort, trait, covset, r["h2"], r["se"], r["n"]))
  res[[length(res) + 1]] <<- data.frame(cohort, trait, covset, n = unname(r["n"]),
                                        h2 = unname(r["h2"]), se = unname(r["se"]))
}

## 2020 cohort (n = 807)
cat("=========== 2020 cohort ===========\n")
K <- read_rel("reanalysis/01_structure/grm")
s <- read.delim("reanalysis/00_data/samples.tsv", stringsAsFactors=FALSE) %>%
  slice(match(rownames(K), clean_id))
stopifnot(identical(s$clean_id, rownames(K)))
K20 <- pd(K)
X0  <- model.matrix(~ 1, s)
XR  <- model.matrix(~ region + sex_raw + age, s)
XRS <- model.matrix(~ region + sex_raw + age + score, s)
add("2020", "binary", "none",                              h2_se(s$binary, K20, X0))
add("2020", "binary", "+ sex/age/region",                  h2_se(s$binary, K20, XR))
add("2020", "binary", "+ sex/age/region + tissue-score",   h2_se(s$binary, K20, XRS))
add("2020", "score",  "none",                              h2_se(s$score,  K20, X0))
add("2020", "score",  "+ sex/age/region",                  h2_se(s$score,  K20, XR))

## 2022 cohort (n = 1,170 with no covariates, 1,166 with sex/age/region, 1,151 with sex/age/farm)
## farm models use farms with >= 10 deer
cat("\n=========== 2022 cohort ===========\n")
K   <- read_rel("reanalysis2022/01_structure/grm2022")
cov <- read.csv("reanalysis2022/data/cov_1170.csv", stringsAsFactors=FALSE) %>%
  mutate(cwd = suppressWarnings(as.numeric(binary)), sexn = suppressWarnings(as.numeric(sex)))
big <- cov %>% filter(!is.na(cwd)) %>% count(farm) %>% filter(n >= 10) %>% pull(farm)
fit22 <- function(d, form){
  Ks <- pd(K[d$CMSID, d$CMSID])
  h2_se(d$cwd, Ks, model.matrix(form, d))
}
add("2022", "binary", "none",             fit22(cov %>% filter(!is.na(cwd)), ~ 1))
add("2022", "binary", "+ sex/age/region", fit22(cov %>% filter(!is.na(cwd), !is.na(sexn), !is.na(age)),
                                                ~ factor(sexn) + age + region))
add("2022", "binary", "+ sex/age/farm",   fit22(cov %>% filter(!is.na(cwd), !is.na(sexn), !is.na(age), farm %in% big),
                                                ~ factor(sexn) + age + factor(farm)))

## write table
out <- bind_rows(res)
write.csv(out, file.path(OUT, "h2_rrblup_se.csv"), row.names = FALSE)
cat("\nwrote", file.path(OUT, "h2_rrblup_se.csv"), "\n")
