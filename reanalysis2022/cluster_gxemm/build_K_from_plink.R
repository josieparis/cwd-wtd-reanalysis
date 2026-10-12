# Build the GxEMM kinship files (K_*.bin, float64) from the PLINK 1.9 --make-rel GRMs, the GRM used in the rest of the paper
# 2022 (n = 1171): grm2022.rel to cluster_gxemm/K_1170.bin; 2020 (n = 807): grm.rel to reanalysis/01_structure/K_807.bin
# the previous NumPy bin is kept as *.numpy.bak

## working dir
setwd("<dir>/Seabury_data")

read_plink_rel <- function(prefix){                    # prefix without extension
  relf <- if (file.exists(paste0(prefix, ".rel.gz"))) gzfile(paste0(prefix, ".rel.gz")) else paste0(prefix, ".rel")
  K   <- as.matrix(read.table(relf, sep = "\t"))
  ids <- read.table(paste0(prefix, ".rel.id"), stringsAsFactors = FALSE)$V2
  stopifnot(nrow(K) == ncol(K), nrow(K) == length(ids))
  dimnames(K) <- list(ids, ids)
  K
}

## align = "id": reorder the GRM to the covariate file by CMSID (2022)
## align = "position": GRM rows already match the covariate rows (2020, IDs renamed S0001.. in the bed),
## checked against the previous bin to catch a misorder
build_one <- function(rel_prefix, cov_csv, out_bin, n_expect, align = "id"){
  cat(sprintf("\n=== %s ===\n", basename(out_bin)))
  K   <- read_plink_rel(rel_prefix)
  cov <- read.csv(cov_csv, stringsAsFactors = FALSE)
  id_col <- if ("CMSID" %in% names(cov)) "CMSID" else names(cov)[1]
  ids_cov <- as.character(cov[[id_col]])
  stopifnot(length(ids_cov) == n_expect, nrow(K) == n_expect)
  if (align == "id"){
    if (!all(ids_cov %in% rownames(K))) stop("cov IDs missing from PLINK GRM: ",
        paste(head(setdiff(ids_cov, rownames(K))), collapse = ", "))
    Kp <- K[ids_cov, ids_cov]                          # covariate file order
    stopifnot(isTRUE(all.equal(ids_cov, rownames(Kp))))
  } else {
    Kp <- unname(K)                                     # positional: row i == covariate row i
  }
  stopifnot(nrow(Kp) == n_expect)

  ## compare with the previous bin, if there is one
  if (file.exists(out_bin)){
    old <- matrix(readBin(out_bin, "double", n_expect * n_expect), n_expect, n_expect, byrow = TRUE)
    lt  <- lower.tri(Kp)
    r   <- cor(Kp[lt], old[lt]); md <- max(abs(Kp[lt] - old[lt]))
    cat(sprintf("  concordance vs previous bin: r = %.6f   max|diff| = %.4f   (diag: PLINK %.3f / prev %.3f)\n",
                r, md, mean(diag(Kp)), mean(diag(old))))
    if (align == "position" && r < 0.999)
      stop("positional alignment failed concordance guard (r < 0.999) — 2020 GRM/cov order may differ")
    bak <- sub("\\.bin$", ".numpy.bak", out_bin)
    if (!file.exists(bak)) { file.copy(out_bin, bak); cat("  backed up previous bin -> ", basename(bak), "\n") }
  }
  writeBin(as.double(t(Kp)), out_bin)                  # row-major float64
  cat(sprintf("  wrote %s  (%d x %d, PLINK --make-rel VanRaden)\n", out_bin, n_expect, n_expect))
  invisible(Kp)
}

## build the 2022 and 2020 bins
build_one("reanalysis2022/01_structure/grm2022", "reanalysis2022/data/cov_1170.csv",
          "reanalysis2022/cluster_gxemm/K_1170.bin", 1171)
build_one("reanalysis/01_structure/grm",          "reanalysis/01_structure/cov_807.csv",
          "reanalysis/01_structure/K_807.bin",       807, align = "position")

## copy the 2020 inputs into cluster_gxemm_2020 (run_gxemm_2020.R reads them from there)
file.copy("reanalysis/01_structure/K_807.bin",  "reanalysis2022/cluster_gxemm_2020/K_807.bin",  overwrite = TRUE)
file.copy("reanalysis/01_structure/cov_807.csv", "reanalysis2022/cluster_gxemm_2020/cov_807.csv", overwrite = TRUE)
cat("\nstaged K_807.bin + cov_807.csv into cluster_gxemm_2020/\nDone.\n")
