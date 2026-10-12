# ARSB intron 4 SNP (Affx-573985690): MAF and CWD prevalence per farm for the 16 farms with >= 10 deer
# (Supplementary Fig 3), and the correlation across farms within each region (Results section 5)

library(dplyr)

## working dir
setwd("<dir>/Seabury_data")

## read in covariates and candidate SNP genotypes
cv <- read.csv("reanalysis2022/data/cov_1170.csv", stringsAsFactors = FALSE)
gt <- read.delim("reanalysis2022/data/candidate_genos_1170.tsv", check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(gt$CMSID == cv$CMSID))
cv <- cv %>% mutate(y = suppressWarnings(as.numeric(binary)), g = suppressWarnings(as.numeric(gt[["Affx-573985690"]])))
p <- mean(cv$g, na.rm = TRUE) / 2; if (p > 0.5) cv <- cv %>% mutate(g = 2 - g)   # count the minor allele

## prevalence and MAF per farm (farms with >= 10 deer of known status)
out <- cv %>% filter(!is.na(y)) %>% group_by(farm) %>%
  summarise(region = first(region), n = n(), prev = mean(y), maf = mean(g, na.rm = TRUE) / 2, .groups = "drop") %>%
  filter(n >= 10) %>% arrange(region, farm) %>% as.data.frame()

## write table
write.csv(out, "reanalysis2022/data/arsb_perfarm.csv", row.names = FALSE)
print(out, row.names = FALSE)

## correlation of MAF and prevalence across farms, per region
for (r in unique(out$region)) { d <- out %>% filter(region == r); ct <- cor.test(d$maf, d$prev)
  cat(sprintf("%-8s k = %d farms: r = %+.2f (P = %.3f)\n", r, nrow(d), ct$estimate, ct$p.value)) }
