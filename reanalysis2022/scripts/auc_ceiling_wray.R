# AUC ceiling (Wray et al. 2010) and the fraction of genetic variance (R) a predictor needs to reach each observed
# AUC, at prevalence K = 0.01 to 0.35 using the liability h2 at each K (Lee et al. 2011 transform)
# run after prediction_lofo.R, writes MANUSCRIPT_tables/auc_ceiling_wray.csv

## working dir
setwd("<dir>/Seabury_data")

## AUC of a genetic predictor that captures a fraction R of the liability genetic variance (R = 1 is the ceiling)
auc_liab <- function(K, h2L, R = 1) {                # Wray et al. 2010
  V  <- R * h2L; Th <- qnorm(1 - K); z <- dnorm(Th)
  i  <- z / K; ip <- -z / (1 - K)                    # mean liability of cases / controls
  vU <- 1 - i * (i - Th); vL <- 1 - ip * (ip - Th)   # variance reduction among cases / controls
  pnorm((V * i - V * ip) / sqrt(V^2 * vU + V * (1 - V) + V^2 * vL + V * (1 - V)))
}
## R needed to reach an observed AUC
implied_R <- function(K, h2L, auc) {                 # NA if above the ceiling
  if (h2L >= 1 || auc >= auc_liab(K, h2L, 1)) return(NA_real_)
  uniroot(function(R) auc_liab(K, h2L, R) - auc, c(1e-8, 1))$root
}
## observed to liability scale h2 (Lee et al. 2011)
lee <- function(h2o, P, K) { z <- dnorm(qnorm(1 - K)); h2o * (K * (1 - K))^2 / (P * (1 - P) * z^2) }
## read a value from a GCTA .hsq file
hsq <- function(f, key = "V(G)/Vp") { L <- strsplit(readLines(f), "\t"); as.numeric(Filter(function(z) z[1] == key, L)[[1]][2]) }

## check against Tsairidou et al. 2014 (K = 0.5, h2_L = 0.34: AUC_max 0.77, AUC 0.69 at R = 0.5)
stopifnot(abs(auc_liab(0.5, 0.34) - 0.77) < 0.005, abs(auc_liab(0.5, 0.34, 0.5) - 0.69) < 0.005)

## observed scale h2 (GCTA, PLINK GRM) and case fraction of each sample
G22 <- "reanalysis2022/cluster_h2/gcta_plink/results"; G20 <- "reanalysis2022/cluster_h2/gcta_2020/results"
cv  <- read.csv("reanalysis2022/data/cov_1170.csv", stringsAsFactors = FALSE)
y   <- suppressWarnings(as.numeric(cv$binary)); sx <- suppressWarnings(as.numeric(cv$sex))
tb  <- table(cv$farm[!is.na(y)]); sF <- !is.na(y) & !is.na(sx) & !is.na(cv$age) & cv$farm %in% names(tb)[tb >= 10]
mods <- list(
  `2022, unadjusted`    = list(h2o = hsq(file.path(G22, "h2_nocov.hsq")), P = mean(y, na.rm = TRUE)),
  `2022, farm-adjusted` = list(h2o = hsq(file.path(G22, "h2_farm.hsq")),  P = mean(y[sF])),
  `2020, unadjusted`    = list(h2o = hsq(file.path(G20, "b_none.hsq")),   P = 284 / 807))
# check the transform against GCTA --prevalence 0.05
stopifnot(abs(lee(mods[[1]]$h2o, mods[[1]]$P, 0.05) - hsq(file.path(G22, "h2_nocov_liab_0.05.hsq"), "V(G)/Vp_L")) < 0.002,
          abs(lee(mods[[3]]$h2o, mods[[3]]$P, 0.05) - hsq(file.path(G20, "b_none_liab_0.05.hsq"), "V(G)/Vp_L")) < 0.002)

## observed AUCs: 2022 random 5-fold and leave-one-farm-out (Table 2), 2020 in-sample (Supp Table 4)
pr  <- readRDS("reanalysis2022/data/prediction_lofo_plink.rds")$none
t20 <- read.csv("reanalysis/01_structure/table1_reproduction.csv", stringsAsFactors = FALSE)
AUC <- c(in22 = mean(pr$random_5fold[, "auc"]), lofo22 = mean(pr$lofo[, "auc"]),
         in20 = t20$mean[t20$config == "k=5  None" & t20$metric == "AUC"])

## liability h2, AUC ceiling and implied R at each K
KS  <- c(0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.26, 0.30, 0.35)
out <- do.call(rbind, lapply(names(mods), function(m) do.call(rbind, lapply(KS, function(K) {
  h <- lee(mods[[m]]$h2o, mods[[m]]$P, K); a <- if (grepl("^2020", m)) AUC["in20"] else AUC[c("in22", "lofo22")]
  data.frame(model = m, K = K, h2L = round(h, 3), AUC_max = round(if (h < 1) auc_liab(K, h) else NA, 3),
             R_insample = round(implied_R(K, h, a[1]), 3),
             R_out_of_farm = if (length(a) > 1) round(implied_R(K, h, a[2]), 4) else NA)
}))))

## write table
dir.create("MANUSCRIPT_tables", showWarnings = FALSE)
write.csv(out, "MANUSCRIPT_tables/auc_ceiling_wray.csv", row.names = FALSE)

## summary
cat(sprintf("observed AUC: 2022 in-sample %.3f | 2022 out-of-farm %.3f | 2020 in-sample %.3f\n", AUC[1], AUC[2], AUC[3]))
for (m in names(mods)) cat(sprintf("%-20s h2_obs = %.3f, case fraction = %.3f\n", m, mods[[m]]$h2o, mods[[m]]$P))
print(out, row.names = FALSE)
rng <- function(x) sprintf("%.2f to %.2f", min(x, na.rm = TRUE), max(x, na.rm = TRUE))
for (lo in c(0.01, 0.05)) {
  s <- out[out$K >= lo, ]
  cat(sprintf("\nK = %.2f to 0.35:\n", lo))
  for (m in names(mods)) { d <- s[s$model == m, ]
    cat(sprintf("  %-20s AUC_max %s | in-sample needs R %s (%d K above the ceiling)%s\n", m, rng(d$AUC_max), rng(d$R_insample),
                sum(is.na(d$R_insample)), if (!grepl("^2020", m)) sprintf(" | out-of-farm needs R %s", rng(d$R_out_of_farm)) else "")) }
}
cat("\nwrote MANUSCRIPT_tables/auc_ceiling_wray.csv\n")
