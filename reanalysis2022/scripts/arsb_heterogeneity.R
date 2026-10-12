# Results section 5: ARSB intron-4 SNP in the authors' METAL meta-analyses (Dryad files AF22, AF25, AF28)
# Cochran's Q heterogeneity P for ARSB, and how many SNPs pass the authors' polygenic threshold (P < 5e-5)
# compared with the number expected by chance

library(readxl)
library(dplyr)

## working dir
setwd("<dir>/Seabury_data")

ZIP  <- "doi_10_5061_dryad_wh70rxwnt__v20220622.zip"
W    <- "reanalysis2022/work"
ARSB <- "Affx-573985690"
QP   <- "P-Value from Cochran's Q"
files <- c(none = "AF22_", sex_age = "AF25_", sex_age_farm = "AF28_")

for (m in names(files)) {
  f <- grep(files[[m]], unzip(ZIP, list = TRUE)$Name, value = TRUE)
  if (!file.exists(file.path(W, f))) unzip(ZIP, files = f, exdir = W)
  x <- suppressWarnings(read_excel(file.path(W, f), .name_repair = "unique_quiet")) %>%
    filter(`Studies for which Valid` == 3) %>%                  # SNPs tested in all three regions
    transmute(marker = Marker, qp = suppressWarnings(as.numeric(.data[[QP]])),
              main = suppressWarnings(as.numeric(`Overall P-Value (Main Effect)`))) %>%
    filter(!is.na(qp))
  a <- filter(x, marker == ARSB)
  lam <- median(qchisq(x$qp, 2, lower.tail = FALSE)) / qchisq(0.5, 2)
  cat(sprintf("%-13s SNPs %s | lambda(Q) %.2f | Q P < 5e-5: %d (expected by chance %.1f) | ARSB: Q P = %.2g (rank %d), main-effect P = %.2f\n",
              m, format(nrow(x), big.mark = ","), lam, sum(x$qp < 5e-5), nrow(x) * 5e-5,
              a$qp, sum(x$qp <= a$qp), a$main))
}
