# helper sourced by the 2022 scripts: load_grm2022() reads the PLINK GRM (plink_2022.sh) in the order of ids
# which = "all" (119,962 SNPs) or "no13" (without chr13, which carries PRNP; 115,272 SNPs)
# bend = TRUE floors the eigenvalues at 1e-6 so the GRM is positive definite for rrBLUP
load_grm2022 <- function(ids, which = c("all", "no13"), bend = TRUE) {
  which <- match.arg(which)
  pre <- c(all = "reanalysis2022/01_structure/grm2022", no13 = "reanalysis2022/01_structure/grm2022_no13")[[which]]
  K   <- as.matrix(read.table(gzfile(paste0(pre, ".rel.gz"))))
  kid <- read.table(paste0(pre, ".rel.id"), stringsAsFactors = FALSE)$V2
  dimnames(K) <- list(kid, kid)
  stopifnot(all(ids %in% kid))
  K <- K[ids, ids]
  if (bend) {
    e <- eigen(K, symmetric = TRUE)
    K <- e$vectors %*% (pmax(e$values, 1e-6) * t(e$vectors))
    dimnames(K) <- list(ids, ids)
  }
  K
}
