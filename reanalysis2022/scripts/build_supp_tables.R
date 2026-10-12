# Supplementary Tables 1-4 (Supplementary Table 5, per-farm GxEMM heritability, is made by build_table5.R)
# S1 APHIS-VS funding for CWD predictive genetics 2021-2025, S2 the 2020 and 2022 datasets,
# S3 naive and 10-PC association scans, S4 reproduction of the GBLUP accuracy in Seabury et al. (2020) Table 1

library(dplyr)

## working dir
setwd("<dir>/Seabury_data")

OUT <- "MANUSCRIPT_tables"

## formatting helpers
f0  <- function(x) formatC(x, format = "d", big.mark = ",")
usd <- function(x) paste0("$", formatC(x, format = "f", digits = 2, big.mark = ","))
sup <- c("0"="⁰","1"="¹","2"="²","3"="³","4"="⁴","5"="⁵","6"="⁶","7"="⁷","8"="⁸","9"="⁹","-"="⁻")
sci <- function(p){ e <- floor(log10(p))
  sprintf("%.1f×10%s", p / 10^e, paste(sup[strsplit(as.character(e), "")[[1]]], collapse = "")) }
## genomic control lambda
lam <- function(p){ p <- p[is.finite(p) & p > 0]; median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1) }
lamf <- function(l) if (l >= 2) sprintf("%.1f", l) else sprintf("%.2f", l)
wr  <- function(d, f){ write.csv(d, file.path(OUT, f), row.names = FALSE); cat("wrote", file.path(OUT, f), "\n") }

## S1: APHIS-VS predictive genetics funding
## transcribed from the APHIS-VS table "CA Funds Awarded to Predictive Genetics" (2025 awarded non-competitively)
## the source's Total Funding cell says $22,528,638.20 but the yearly rows sum to $22,582,638.20, used here
s1 <- data.frame(Year = c("2021", "2022", "2023", "2024", "2025*"),
                 pg   = c(1228946.68, 1916848.31, 1785803, 3006338, 2796862),
                 tot  = c(2775000, 4755668.20, 5258473, 6045215, 3748282),
                 org  = c("State Ag (4)", "State Ag (4)", "State Ag (8); Univ (1)",
                          "State Ag (15); Univ (1)", "State Ag (16); Univ (1)"))
stopifnot(isTRUE(all.equal(sum(s1$pg), 10734797.99)))           # the source's total
S1 <- bind_rows(
  s1 %>% transmute(Year, `Predictive-genetics funding` = usd(pg),
                   `Total funding` = usd(tot), `Percent of total` = sprintf("%.0f%%", 100 * pg / tot),
                   `Organizations awarded` = org),
  s1 %>% summarise(Year = "Total", `Predictive-genetics funding` = usd(sum(pg)),
                   `Total funding` = usd(sum(tot)), `Percent of total` = sprintf("%.0f%%", 100 * sum(pg) / sum(tot)),
                   `Organizations awarded` = "—"))
wr(S1, "tableS1_aphis_funding.csv")

## S2: the two datasets
s20 <- read.delim("reanalysis/00_data/samples.tsv", stringsAsFactors = FALSE)
c22 <- read.csv("reanalysis2022/data/cov_1170.csv", stringsAsFactors = FALSE) %>%
  mutate(cwd = suppressWarnings(as.numeric(binary)), sexn = suppressWarnings(as.numeric(sex)))
k22 <- c22 %>% filter(!is.na(cwd))
relab <- function(r) if_else(r == "East", "Northeast", r)
REG <- c("Northeast", "Midwest", "South")
rn  <- function(reg) paste(as.integer(table(factor(relab(reg), REG))), collapse = " / ")
rp  <- function(reg, y) paste(sprintf("%.2f", tapply(y, factor(relab(reg), REG), mean)), collapse = " / ")
fm22 <- k22 %>% group_by(farm) %>% summarise(n = n(), both = n_distinct(cwd) > 1)   # deer per farm; both = has cases and controls
big  <- fm22 %>% filter(n >= 10) %>% pull(farm)
cc22 <- k22 %>% filter(!is.na(sexn), !is.na(age))
n_f10  <- sum(cc22$farm %in% big)
n_f10b <- sum(cc22$farm %in% (fm22 %>% filter(n >= 10, both) %>% pull(farm)))
nbim   <- function(f) length(readLines(f))
auto20 <- sum(read.table("reanalysis/00_data/deer.bim")$V1 %in% 1:29)          # GRM set (grm.log)
assoc20 <- nbim("reanalysis/02_assoc/deer_bin.bim")                            # MAF >= 0.05 association set
grm22  <- nbim("reanalysis2022/01_structure/deer2022_qc.bim")                  # after --geno 0.05
assoc22 <- nrow(read.delim("reanalysis2022/mlma_results/mlma_bin.mlma"))      # GEMMA / GCTA scans
dups <- function(f){ K <- as.matrix(read.table(gzfile(f))); sum(K[upper.tri(K)] >= 0.9) }
mx  <- read.csv("reanalysis2022/data/figS1_max_per_animal.csv")$maxconc       # best 2022 match per 2020 deer
S2 <- data.frame(check.names = FALSE,
  ` ` = c("Publication", "Dryad deposit", "Deer with CWD status", "Binary CWD (cases / controls)",
          "Sample CWD prevalence", "Ordinal tissue score (0/1/2)", "Covariates deposited", "Farms",
          "Deer by region (Northeast / Midwest / South)", "CWD prevalence by region (Northeast / Midwest / South)",
          "Analysis samples", "SNPs (as reported by the original study)",
          "Autosomal SNPs: relationship matrix (PLINK --make-rel)", "Autosomal SNPs: association scans",
          "Genotype-identical pairs within the dataset", "Overlap with the other dataset"),
  `2020 dataset` = c("Seabury et al. 2020, G3 10:1433", "doi:10.5061/dryad.xd2547dcw", f0(nrow(s20)),
          paste(sum(s20$binary == 1), "/", sum(s20$binary == 0)), sprintf("%.3f", mean(s20$binary)), "yes",
          "sex, age, region", "not deposited", rn(s20$region), rp(s20$region, s20$binary),
          sprintf("%s (all analyses)", f0(nrow(s20))), "123,987", f0(auto20),
          sprintf("%s (minor-allele frequency ≥ 0.05)", f0(assoc20)),
          as.character(dups("reanalysis/01_structure/grm.rel.gz")), "—"),
  `2022 dataset` = c("Seabury et al. 2022, G3 12:jkac109", "doi:10.5061/dryad.wh70rxwnt",
          sprintf("%s (%s deposited)", f0(nrow(k22)), f0(nrow(c22))),
          paste(sum(k22$cwd == 1), "/", sum(k22$cwd == 0)), sprintf("%.3f", mean(k22$cwd)), "no",
          "sex, age, region, farm", sprintf("%d (%d with ≥ 10 deer)", nrow(fm22), length(big)),
          rn(k22$region), rp(k22$region, k22$cwd),
          sprintf("%s (no covariates); %s (complete sex and age); %s (farm models: farms with ≥ 10 deer); %s (GxEMM and LDAK-PCGC: farms with ≥ 10 deer and both cases and controls)",
                  f0(nrow(k22)), f0(nrow(cc22)), f0(n_f10), f0(n_f10b)),
          "125,590", f0(grm22), f0(assoc22), as.character(dups("reanalysis2022/01_structure/grm2022.rel.gz")),
          sprintf("%d of the %d 2020 deer (%.0f%%) are genotype-identical (IBS ≥ 0.90) to a 2022 deer",
                  sum(mx >= 0.90), length(mx), 100 * mean(mx >= 0.90))))
wr(S2, "tableS2_cohorts.csv")

## S3: naive and PC-adjusted scans
PRNP <- "G96S_G_to_A_at286"
g20b <- read.delim(gzfile("reanalysis/01_structure/gwas_binary.tsv.gz"))     # 2020, 120,113 SNPs (Table 1 set)
g20s <- read.delim(gzfile("reanalysis/01_structure/gwas_score.tsv.gz"))
g22  <- read.delim(gzfile("reanalysis2022/data/gwas2022_binary.tsv.gz"))     # 2022, 119,962 SNPs (GRM QC set)
row3 <- function(ds, tr, mod, d, col) { p <- d[[col]]
  c(ds, tr, mod, lamf(lam(p)), f0(sum(p < 5e-8, na.rm = TRUE)), f0(sum(p < 5e-5, na.rm = TRUE)), sci(p[d$marker == PRNP])) }
S3 <- as.data.frame(rbind(
  row3(sprintf("2020 (n = %s)", f0(nrow(s20))), "Binary", "Structure-naive", g20b, "p_naive"),
  row3("", "", "+ 10 genomic PCs", g20b, "p_pc"),
  row3("", "Score", "Structure-naive", g20s, "p_naive"),
  row3("", "", "+ 10 genomic PCs", g20s, "p_pc"),
  row3(sprintf("2022 (n = %s)", f0(nrow(k22))), "Binary", "Structure-naive", g22, "p_naive"),
  row3("", "", "+ 10 genomic PCs", g22, "p_pc")))
names(S3) <- c("Dataset", "Phenotype", "Model", "λ", "SNPs P < 5×10⁻⁸", "SNPs P < 5×10⁻⁵", "PRNP P")
wr(S3, "tableS3_naive_pc.csv")

## S4: reproduction of Seabury et al. (2020) Table 1
rep <- read.csv("reanalysis/01_structure/table1_reproduction.csv", stringsAsFactors = FALSE)   # from 12_table1_reproduce.R
MET <- c("AUC", "Accuracy", "Sensitivity", "Specificity", "MCC", "RMSE")
orig <- rbind(                                               # as reported in Seabury et al. 2020 Table 1
  "k=3  None"           = c(AUC = .8471, Accuracy = .8167, Sensitivity = .6447, Specificity = .9101, MCC = .5871, RMSE = .3768),
  "k=3  Sex,Age,Region" = c(AUC = .8534, Accuracy = .8119, Sensitivity = .6735, Specificity = .8870, MCC = .5787, RMSE = .3746),
  "k=5  None"           = c(AUC = .8485, Accuracy = .8198, Sensitivity = .6496, Specificity = .9121, MCC = .5940, RMSE = .3754),
  "k=5  Sex,Age,Region" = c(AUC = .8542, Accuracy = .8159, Sensitivity = .6716, Specificity = .8942, MCC = .5870, RMSE = .3734))
S4 <- NULL; dmax <- 0
for (cf in rownames(orig)) {
  r <- rep %>% filter(config == cf); m <- setNames(r$mean, r$metric)[MET]; s <- setNames(r$sd, r$metric)[MET]
  dmax <- max(dmax, abs(m - orig[cf, MET]))
  k <- sub("k=([0-9]).*", "\\1", cf); cv <- if (grepl("None", cf)) "none" else "sex, age, region"
  S4 <- rbind(S4,
    c(k, cv, "This study (mean ± SD, 50 replicates)", sprintf("%.3f ± %.3f", m, s)),
    c("", "", "Seabury et al. 2020 (as reported)", sprintf("%.4f", orig[cf, MET])))
}
S4 <- as.data.frame(S4); names(S4) <- c("Folds (k)", "Covariates", "Source", MET)
wr(S4, "tableS4_prediction_repro.csv")
cat(sprintf("S4: max |this study - original| over 24 cells = %.4f\n", dmax))
