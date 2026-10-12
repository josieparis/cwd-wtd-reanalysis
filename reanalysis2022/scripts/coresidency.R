# statistics for Results section 1 (CWD status is confounded with region, farm and relatedness) in both cohorts,
# using the PLINK GRMs and PCs, with PC1 flipped so South is positive as in Fig 1a

library(dplyr)

## working dir
setwd("<dir>/Seabury_data")

## helper functions
sq     <- function(f) as.matrix(read.table(gzfile(f), sep = "\t"))
flipS  <- function(x, reg) if (mean(x[reg == "South"]) < 0) -x else x
relabR <- function(r) ifelse(r == "East", "Northeast", r)
nagel  <- function(m, m0) {                                  # Nagelkerke R2
  n <- length(m$y); l0 <- as.numeric(logLik(m0)); l1 <- as.numeric(logLik(m))
  (1 - exp((2 / n) * (l0 - l1))) / (1 - exp((2 / n) * l0))
}
orp <- function(m, term) { s <- summary(m)$coefficients
  sprintf("OR = %.2f (95%% CI %.2f-%.2f), P = %.2g", exp(s[term, 1]),
          exp(s[term, 1] - 1.96 * s[term, 2]), exp(s[term, 1] + 1.96 * s[term, 2]), s[term, 4]) }

## PC1 vs region and CWD status, and logistic models of CWD on region and 10 PCs
structure_block <- function(tag, pcs, region, y, ev) {
  cat(sprintf("\n================ %s ================\n", tag))
  ok <- !is.na(y) & !is.na(region)
  pc1 <- flipS(pcs[, 1], region)
  cat(sprintf("n with status = %d | prevalence by region: %s\n", sum(ok),
      paste(sprintf("%s %.0f%%", names(tapply(y[ok], region[ok], mean)),
                    100 * tapply(y[ok], region[ok], mean)), collapse = ", ")))
  cat(sprintf("PC1 = %.1f%% of variance across the top 20 PCs (%.1f%% of total GRM trace not computed here)\n",
      100 * ev[1] / sum(ev[1:20]), NA))
  a  <- summary(lm(pc1[ok] ~ region[ok]))$r.squared
  kw <- kruskal.test(pc1[ok] ~ factor(region[ok]))
  k  <- length(unique(region[ok])); nn <- sum(ok)
  cat(sprintf("PC1 ~ region: eta2(ANOVA) = %.3f | eta2(KW-H) = %.3f | Kruskal-Wallis H = %.1f, P = %.2g\n",
      a, (kw$statistic - k + 1) / (nn - k), kw$statistic, kw$p.value))
  ct <- cor.test(pc1[ok], y[ok])
  cat(sprintf("PC1 ~ CWD status: point-biserial r = %.3f (r2 = %.3f), P = %.2g\n",
      ct$estimate, ct$estimate^2, ct$p.value))
  d  <- data.frame(y = y[ok], region = relevel(factor(region[ok]), ref = "South"), pcs[ok, 1:10])
  pcn <- paste0("PC", 1:10); names(d)[3:12] <- pcn
  m0 <- glm(y ~ 1, binomial, d); mR <- glm(y ~ region, binomial, d)
  mP <- glm(reformulate(pcn, "y"), binomial, d); mRP <- glm(reformulate(c("region", pcn), "y"), binomial, d)
  cat(sprintf("Logistic, region only        : Nagelkerke R2 = %.3f | Midwest vs South %s | Northeast vs South %s\n",
      nagel(mR, m0), orp(mR, "regionMidwest"), orp(mR, "regionNortheast")))
  cat(sprintf("Logistic, 10 PCs only        : Nagelkerke R2 = %.3f\n", nagel(mP, m0)))
  cat(sprintf("Logistic, region + 10 PCs    : Nagelkerke R2 = %.3f | Midwest vs South %s\n",
      nagel(mRP, m0), orp(mRP, "regionMidwest")))
  cat(sprintf("LRT region+PCs vs region     : P = %.2g ;  region+PCs vs PCs : P = %.2g\n",
      anova(mR, mRP, test = "LRT")[2, 5], anova(mP, mRP, test = "LRT")[2, 5]))
}

## relatedness: duplicates, first-degree and close relatives, and whether they share a region or farm
relatedness_block <- function(G, region, y, farm = NULL) {
  n <- nrow(G); lt <- which(lower.tri(G), arr.ind = TRUE); g <- G[lt]
  i <- lt[, 1]; j <- lt[, 2]
  dup <- g >= 0.9; fd <- g >= 0.4; cl <- g >= 0.2
  Gz <- G; diag(Gz) <- 0
  cat(sprintf("GRM n = %d | pairs = %d | duplicate pairs (>=0.9) = %d | first-degree pairs (>=0.4, incl. dups) = %d (%d excl. dups)\n",
      n, length(g), sum(dup), sum(fd), sum(fd & !dup)))
  cat(sprintf("animals with >=1 first-degree relative (>=0.4) = %.1f%% | with >=1 close relative (>=0.2) = %.1f%%\n",
      100 * mean(rowSums(Gz >= 0.4) >= 1), 100 * mean(rowSums(Gz >= 0.2) >= 1)))
  sameR <- region[i] == region[j]
  cat(sprintf("first-degree pairs within the same region = %d / %d (%.1f%%)\n",
      sum(fd & sameR, na.rm = TRUE), sum(fd), 100 * mean(sameR[fd], na.rm = TRUE)))
  both <- y[i] == 1 & y[j] == 1; okp <- !is.na(both) & sameR & !dup
  for (u in c(0.05, 0.2)) {
    rel <- okp & g >= 0.2; unr <- okp & g < u
    cat(sprintf("within-region, both CWD+ : related (0.2<=GRM<0.9) %.1f%% of %d | unrelated (GRM<%.2f) %.1f%% of %d\n",
        100 * mean(both[rel]), sum(rel), u, 100 * mean(both[unr]), sum(unr)))
  }
  if (!is.null(farm)) {
    fk <- !is.na(farm) & farm != ""
    sel <- fd & fk[i] & fk[j]; sf <- farm[i][sel] == farm[j][sel]
    xi <- c(i[sel][!sf], j[sel][!sf])
    cat(sprintf("first-degree pairs with BOTH farms recorded = %d of %d (%d involve an animal with no farm label)\n",
        sum(sel), sum(fd), sum(fd) - sum(sel)))
    cat(sprintf("  same farm = %d (%.1f%%) | different farm = %d pairs, involving %d distinct animals\n",
        sum(sf), 100 * mean(sf), sum(!sf), length(unique(xi))))
    di <- which(dup)
    cat("duplicate pairs (GRM >= 0.9):\n")
    for (k in di) cat(sprintf("  GRM = %.3f | farms %s / %s | status %s / %s\n", g[k],
        farm[i[k]], farm[j[k]], y[i[k]], y[j[k]]))
  }
}

## 2020 cohort
p20 <- read.delim("reanalysis/01_structure/pca_merged.tsv", check.names = FALSE, stringsAsFactors = FALSE)
ev20 <- scan("reanalysis/01_structure/pca.eigenval", quiet = TRUE)
structure_block("2020 cohort", as.matrix(p20[, paste0("PC", 1:10)]), relabR(p20$region),
                suppressWarnings(as.numeric(p20$binary)), ev20)
c807 <- read.csv("reanalysis/01_structure/cov_807.csv", stringsAsFactors = FALSE)
relatedness_block(sq("reanalysis/01_structure/grm.rel.gz"), relabR(c807$region),
                  suppressWarnings(as.numeric(c807$binary)))           # GRM rows in cov_807 order

## 2022 cohort
pe   <- read.table("reanalysis2022/01_structure/pca2022.eigenvec", stringsAsFactors = FALSE)
ev22 <- scan("reanalysis2022/01_structure/pca2022.eigenval", quiet = TRUE)
cov  <- read.csv("reanalysis2022/data/cov_1170.csv", stringsAsFactors = FALSE)
mi   <- match(pe$V2, cov$CMSID)
structure_block("2022 cohort", as.matrix(pe[, 3:12]), relabR(cov$region[mi]),
                suppressWarnings(as.numeric(cov$binary[mi])), ev22)
ids  <- read.table("reanalysis2022/01_structure/grm2022.rel.id", stringsAsFactors = FALSE)$V2
mg   <- match(ids, cov$CMSID)
relatedness_block(sq("reanalysis2022/01_structure/grm2022.rel.gz"), relabR(cov$region[mg]),
                  suppressWarnings(as.numeric(cov$binary[mg])), cov$farm[mg])

## farms with >= 10 deer and their range of prevalence
fp <- cov %>% mutate(y = suppressWarnings(as.numeric(binary))) %>% filter(!is.na(y), !is.na(farm), farm != "") %>%
  group_by(farm) %>% summarise(n = n(), prev = mean(y), .groups = "drop") %>% filter(n >= 10)
cat(sprintf("farms with >=10 deer = %d | prevalence range %.0f%%-%.0f%%\n", nrow(fp), 100*min(fp$prev), 100*max(fp$prev)))
