# Supplementary Figure 1: genotype identity (IBS) between the 2020 and 2022 datasets
# data from figS1_concordance_compute.py: IBS allele sharing for all 807 x 1171 cross-dataset pairs
# (14,825 markers left from a subsample of 20k shared markers after allele alignment)

library(ggplot2)
library(scales)

## working dir
setwd("<dir>/Seabury_data")

source("reanalysis2022/scripts/fig_theme.R"); TH <- theme_seabury()

## read in data
h  <- read.csv("reanalysis2022/data/figS1_hist_pairwise.csv")
mx <- read.csv("reanalysis2022/data/figS1_max_per_animal.csv")$maxconc
THR <- 0.90
nmatch <- sum(mx >= THR); ntot <- length(mx)

hb <- subset(h, count > 0)
hb$grp <- factor(ifelse(hb$mid >= THR, "same animal (IBS ≥ 0.90)", "distinct animals"),
                 levels = c("distinct animals", "same animal (IBS ≥ 0.90)"))
YMAX <- max(hb$count)
ay   <- 10^(0.66 * log10(YMAX))          # annotation height on the log axis

## plot
p <- ggplot(hb, aes(mid, count, fill = grp)) +
  annotate("rect", xmin = THR, xmax = 1.006, ymin = 0.7, ymax = Inf,
           fill = "#fdd0d0", alpha = 0.45) +
  geom_col(width = 0.005) +
  geom_vline(xintercept = THR, linetype = 2, color = "#b2182b", linewidth = 0.7) +
  annotate("segment", x = 0.893, xend = 0.994, y = ay, yend = nmatch * 1.25,
           color = "#1a4f86", linewidth = 0.4,
           arrow = arrow(length = unit(0.16, "cm"), type = "closed")) +
  annotate("text", x = 0.889, y = ay, hjust = 1, vjust = 0.4, size = 4, color = "#1a4f86",
           label = sprintf("%d / %d 2020 animals (%.0f%%)\nreappear in 2022 at median IBS 1.00",
                           nmatch, ntot, 100 * nmatch / ntot)) +
  annotate("text", x = 0.66, y = ay, hjust = 0.5, vjust = 0.4, size = 4, color = "black",
           label = "distinct animals\nmedian IBS 0.67") +
  scale_fill_manual(values = c("distinct animals" = "grey55",
                               "same animal (IBS ≥ 0.90)" = "#2166AC")) +
  scale_y_log10(labels = label_comma(), expand = expansion(mult = c(0, 0.06))) +
  scale_x_continuous(breaks = seq(0.4, 1.0, 0.1), expand = expansion(mult = 0.015)) +
  labs(x = "between-dataset genotype concordance (IBS allele-sharing proportion)",
       y = "pair count (log scale)",
       title = "Genotype identity between the 2020 and 2022 datasets") +
  TH + theme(legend.position = c(0.72, 0.98), legend.justification = c(1, 1))   # clear of the bars

## save
ggsave("Figures/figS1_concordance.png", p, width = 9.2, height = 6.0, dpi = 175, bg = "white", device=ragg::agg_png)
cat(sprintf("wrote Figures/figS1_concordance.png  (nmatch=%d/%d, markers used=14825)\n",
            nmatch, ntot))
