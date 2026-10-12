# Results section 2: PRNP variance explained in the 2022 deer carried over from 2020 vs the deer new in 2022
# 2022 deer are matched to 2020 deer by genotype concordance (as for Supplementary Fig 1), then the PRNP
# codon 96 R2 is compared between the two groups and the difference in effect is tested (genotype x group)

library(dplyr)

## working dir
setwd("<dir>/Seabury_data")
W <- "reanalysis2022/work"   # deer2022 plink files (plink_2022.sh)

## fingerprint markers: autosomal SNPs in both datasets, call rate >= 98%, MAF >= 0.1, 6,000 drawn at random
b20 <- read.table("reanalysis/00_data/deer.bim")
b22 <- read.table(file.path(W, "deer2022.bim"))
writeLines(b20 %>% filter(V1 %in% 1:29, V2 %in% b22$V2) %>% pull(V2), file.path(W, "fp_common.txt"))
system(paste("plink --bfile reanalysis/00_data/deer --chr-set 29 --extract", file.path(W, "fp_common.txt"),
             "--geno 0.02 --maf 0.1 --thin-count 6000 --seed 1 --recode A --out", file.path(W, "fp20")))
r20 <- read.table(file.path(W, "fp20.raw"), header = TRUE, check.names = FALSE)
mk  <- sub("_[^_]*$", "", names(r20)[-(1:6)])
writeLines(mk, file.path(W, "fp_markers.txt"))
system(paste("plink --bfile", file.path(W, "deer2022"), "--chr-set 29 --extract", file.path(W, "fp_markers.txt"),
             "--recode A --out", file.path(W, "fp22")))
r22 <- read.table(file.path(W, "fp22.raw"), header = TRUE, check.names = FALSE)

## genotype matrices on the same markers; 2022 alleles are coded A/B, so orient them by allele frequency
A <- as.matrix(r20[, -(1:6)]); colnames(A) <- mk
B <- as.matrix(r22[, -(1:6)]); colnames(B) <- sub("_[^_]*$", "", colnames(B)); B <- B[, mk]
p20 <- colMeans(A, na.rm = TRUE) / 2; p22 <- colMeans(B, na.rm = TRUE) / 2
flip <- abs(p20 - p22) > abs(p20 - (1 - p22)); B[, flip] <- 2 - B[, flip]
p22 <- colMeans(B, na.rm = TRUE) / 2
keep <- (p20 < 0.42 | p20 > 0.58) & abs(p20 - p22) < 0.1      # unambiguous orientation only
A <- A[, keep]; B <- B[, keep]
cat(sprintf("fingerprint markers used: %d\n", sum(keep)))

## IBS concordance for every 2022 x 2020 pair: 1 - sum|a - b| / (2 x markers called in both)
ind <- function(M, g) { x <- (M == g) * 1; x[is.na(x)] <- 0; x }
called <- ((!is.na(B)) * 1) %*% t((!is.na(A)) * 1)
D <- 0
for (ga in 0:2) for (gb in 0:2) if (ga != gb) D <- D + abs(ga - gb) * (ind(B, gb) %*% t(ind(A, ga)))
conc <- 1 - D / (2 * pmax(called, 1))
best <- apply(conc, 1, max)
m <- tibble(CMSID = r22$IID, best = best, match2020 = if_else(best >= 0.9, r20$IID[apply(conc, 1, which.max)], NA_character_))
cat(sprintf("2022 records with a 2020 match (IBS >= 0.90): %d of %d (%d distinct 2020 deer)\n",
            sum(!is.na(m$match2020)), nrow(m), n_distinct(na.omit(m$match2020))))

## PRNP codon 96 dosage and CWD status in 2022
d <- read.delim("reanalysis2022/data/candidate_genos_1170.tsv", check.names = FALSE, stringsAsFactors = FALSE) %>%
  transmute(CMSID, y = suppressWarnings(as.numeric(Binary)), G96 = suppressWarnings(as.numeric(G96S_G_to_A_at286))) %>%
  left_join(m, by = "CMSID") %>%
  mutate(group = factor(if_else(is.na(match2020), "new", "carried over"), levels = c("carried over", "new"))) %>%
  filter(!is.na(y), !is.na(G96))

## R2 (squared correlation of CWD status with dosage) and per-allele odds ratio by group
or_ci <- function(y, g) { f <- glm(y ~ g, binomial); sprintf("%.2f (%.2f-%.2f)", exp(coef(f)[2]), exp(confint.default(f)[2, 1]), exp(confint.default(f)[2, 2])) }
res <- d %>% group_by(group) %>%
  summarise(n = n(), prevalence = mean(y), freq_S = mean(G96) / 2, R2 = cor(y, G96)^2, OR = or_ci(y, G96))
print(res)
cat(sprintf("all 2022 deer: n = %d, R2 = %.4f\n", nrow(d), cor(d$y, d$G96)^2))

## does the PRNP effect differ between carried-over and new deer?
fit <- glm(y ~ G96 * group, binomial, d)
cat(sprintf("G96 x group interaction: P = %.4f\n", coef(summary(fit))["G96:groupnew", 4]))
