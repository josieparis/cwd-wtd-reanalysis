# Supplementary Figure 1: genotype identity (IBS) between the 2020 and 2022 datasets
# data from figS1_concordance_compute.py: IBS allele sharing for all 807 x 1171 cross-dataset pairs
# (14,825 markers left from a subsample of 20k shared markers after allele alignment)

library(ggplot2)
library(scales)
library(dplyr)

## working dir
setwd("<dir>/Seabury_data")

source("reanalysis2022/scripts/fig_theme.R"); TH <- theme_seabury()

## read in data
h  <- read.csv("reanalysis2022/data/figS1_hist_pairwise.csv")
mx <- read.csv("reanalysis2022/data/figS1_max_per_animal.csv")$maxconc
THR <- 0.90
nmatch <- sum(mx >= THR); ntot <- length(mx)

## legend labels carry the summary numbers
LAB <- c(distinct = "distinct animals (median IBS 0.67)",
         same     = sprintf("same animal (IBS ≥ 0.90): %d of %d 2020 animals (%.0f%%)\nreappear in 2022 at median IBS 1.00",
                            nmatch, ntot, 100 * nmatch / ntot))
hb <- h %>%
  filter(count > 0) %>%
  mutate(grp = factor(if_else(mid >= THR, LAB[["same"]], LAB[["distinct"]]), levels = unname(LAB)))

## plot
p <- ggplot(hb, aes(mid, count, fill = grp)) +
  annotate("rect", xmin = THR, xmax = 1.006, ymin = 0.7, ymax = Inf,
           fill = "#fdd0d0", alpha = 0.25) +
  geom_col(width = 0.005) +
  geom_vline(xintercept = THR, linetype = 2, color = "#b2182b", linewidth = 0.7) +
  scale_fill_manual(values = setNames(c("grey55", "#2166AC"), unname(LAB))) +
  scale_y_log10(labels = label_comma(), expand = expansion(mult = c(0, 0.06))) +
  scale_x_continuous(breaks = seq(0.4, 1.0, 0.05), expand = expansion(mult = 0.015)) +
  labs(x = "between-dataset genotype concordance (IBS allele-sharing proportion)",
       y = "pair count (log scale)") +
  TH + theme(legend.position = "top", legend.direction = "vertical", legend.justification = "left",
             legend.key.spacing.y = unit(4, "pt"))

## save
ggsave("Figures/figS1_concordance.png", p, width = 9.2, height = 6.0, dpi = 175, bg = "white", device=ragg::agg_png)
cat(sprintf("wrote Figures/figS1_concordance.png  (nmatch=%d/%d, markers used=14825)\n",
            nmatch, ntot))
