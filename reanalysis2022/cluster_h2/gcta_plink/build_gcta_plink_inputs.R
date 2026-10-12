# Build the GCTA inputs for the 2022 heritability runs on the PLINK 1.9 --make-rel GRM (1,171 deer, 119,962 SNPs)
# each model uses the original paper's complete-case sample: 1,170 (no covariates), 1,166 (+ sex/age/region)
# and 1,151 (+ sex/age/farm, farms with >= 10 deer, same deer as GEMMA pheno_farm); writes into cluster_h2/gcta_plink

## working dir
setwd("<dir>/Seabury_data")

OUT <- "reanalysis2022/cluster_h2/gcta_plink"
OLD <- "reanalysis2022/cluster_h2/gcta"
NSNP <- 119962

## GRM: PLINK square .rel to GCTA binary (lower triangle with diagonal, row-wise, float32)
ids <- read.table("reanalysis2022/01_structure/grm2022.rel.id", stringsAsFactors = FALSE)
K   <- as.matrix(read.table(gzfile("reanalysis2022/01_structure/grm2022.rel.gz")))
stopifnot(nrow(K) == ncol(K), nrow(K) == nrow(ids), isSymmetric(unname(K), tol = 1e-6))
lt  <- K[upper.tri(K, diag = TRUE)]          # column-major upper = row-wise lower, as GCTA expects
con <- file(file.path(OUT, "deer2022p.grm.bin"), "wb");   writeBin(lt, con, size = 4); close(con)
con <- file(file.path(OUT, "deer2022p.grm.N.bin"), "wb"); writeBin(rep(NSNP, length(lt)), con, size = 4); close(con)
write.table(ids, file.path(OUT, "deer2022p.grm.id"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

## phenotype and covariates (copied from the earlier run in cluster_h2/gcta)
for (f in c("pheno.txt", "covar_region.txt", "qcovar_age.txt"))
  file.copy(file.path(OLD, f), file.path(OUT, f), overwrite = TRUE)

## farm model: complete sex and age, farms with >= 10 deer of known CWD status
cov <- read.csv("reanalysis2022/data/cov_1170.csv", stringsAsFactors = FALSE)
cov$cwd <- suppressWarnings(as.numeric(cov$binary)); cov$sexn <- suppressWarnings(as.numeric(cov$sex))
tb  <- table(cov$farm[!is.na(cov$cwd)]); big <- names(tb)[tb >= 10]
f10 <- cov[!is.na(cov$cwd) & !is.na(cov$sexn) & !is.na(cov$age) & cov$farm %in% big, ]
g   <- scan("reanalysis2022/cluster_h2/gemma/pheno_farm.txt", what = "", quiet = TRUE)
stopifnot(nrow(f10) == 1151, identical(cov$CMSID[!is.na(g)], f10$CMSID))   # same deer as GEMMA
write.table(data.frame(f10$CMSID, f10$CMSID), file.path(OUT, "keep_farm10.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(data.frame(f10$CMSID, f10$CMSID, f10$sexn, f10$farm), file.path(OUT, "covar_farm10.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

## check the binary GRM reads back to the PLINK matrix
n  <- nrow(K); rb <- readBin(file.path(OUT, "deer2022p.grm.bin"), "numeric", n * (n + 1) / 2, size = 4)
cat(sprintf("GRM %d x %d written (%d values); max |readback - PLINK| = %.2g; farm set n = %d on %d farms\n",
            n, n, length(rb), max(abs(rb - lt)), nrow(f10), length(unique(f10$farm))))
