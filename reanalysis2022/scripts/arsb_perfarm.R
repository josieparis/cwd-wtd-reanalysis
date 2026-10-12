# ARSB intron 4 SNP (Affx-573985690): MAF and CWD prevalence per farm for the 16 farms with >= 10 deer
# (Supplementary Fig 3), and the correlation across farms within each region (Results section 5)

## working dir
setwd("<dir>/Seabury_data")

## read in covariates and candidate SNP genotypes
cv <- read.csv("reanalysis2022/data/cov_1170.csv", stringsAsFactors = FALSE)
gt <- read.delim("reanalysis2022/data/candidate_genos_1170.tsv", check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(gt$CMSID == cv$CMSID))
y <- suppressWarnings(as.numeric(cv$binary)); g <- suppressWarnings(as.numeric(gt[["Affx-573985690"]]))
p <- mean(g, na.rm = TRUE) / 2; if (p > 0.5) g <- 2 - g                     # count the minor allele

## prevalence and MAF per farm
tb <- table(cv$farm[!is.na(y)]); big <- names(tb)[tb >= 10]
out <- do.call(rbind, lapply(big, function(f) { i <- cv$farm == f & !is.na(y)
  data.frame(farm = f, region = cv$region[i][1], n = sum(i), prev = mean(y[i]), maf = mean(g[i], na.rm = TRUE) / 2) }))
out <- out[order(out$region, out$farm), ]

## write table
write.csv(out, "reanalysis2022/data/arsb_perfarm.csv", row.names = FALSE)
print(out, row.names = FALSE)

## correlation of MAF and prevalence across farms, per region
for (r in unique(out$region)) { d <- out[out$region == r, ]; ct <- cor.test(d$maf, d$prev)
  cat(sprintf("%-8s k = %d farms: r = %+.2f (P = %.3f)\n", r, nrow(d), ct$estimate, ct$p.value)) }
