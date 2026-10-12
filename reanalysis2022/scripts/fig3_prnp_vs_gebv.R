# Figure 3: out-of-farm AUC per farm for PRNP codon 96 alone, the full GEBV
# and the GEBV without the PRNP chromosome (polygenic remainder)

library(tidyverse)
library(ggtext)

## working dir
setwd("<dir>/Seabury_data")

source("reanalysis2022/scripts/fig_theme.R"); TH <- theme_seabury()
CPRNP <- "#1B7837"; CGEBV <- "#B2182B"; CREM <- "#E08214"

## read in per-farm AUCs
d <- read.csv("reanalysis2022/data/prnp_threeway_perfarm.csv", stringsAsFactors=FALSE) %>%
  mutate(farm=gsub("East","Northeast", farm)) %>%
  arrange(n) %>%                                                   # order farms by sample size
  mutate(star=pmax(auc_GEBVfull, auc_GEBVminus13) > auc_PRNP,     # asterisk if a GEBV beats PRNP
         lab =sprintf("%s (n=%d)%s", farm, n, if_else(star, "*", "")),
         y   =row_number())
mP <- mean(d$auc_PRNP); mG <- mean(d$auc_GEBVfull); mR <- mean(d$auc_GEBVminus13)
lv <- c(sprintf("*PRNP* codon-96 alone (mean %.2f)", mP),
        sprintf("full genome-wide GEBV (mean %.2f)", mG),
        sprintf("GEBV minus *PRNP* chromosome (mean %.2f)", mR))
dl <- d %>%
  pivot_longer(c(auc_PRNP, auc_GEBVfull, auc_GEBVminus13), names_to="pred", values_to="auc", cols_vary="slowest") %>%
  transmute(y, auc, pred=factor(pred, levels=c("auc_PRNP", "auc_GEBVfull", "auc_GEBVminus13"), labels=lv))
seg <- d %>% transmute(y, xmin=pmin(auc_PRNP,auc_GEBVfull,auc_GEBVminus13),
                       xmax=pmax(auc_PRNP,auc_GEBVfull,auc_GEBVminus13))
cols <- setNames(c(CPRNP,CGEBV,CREM), lv)

## tied AUCs are drawn as one point split between the colours
dl <- dl %>% mutate(key=paste(y, auc), tied=key %in% key[duplicated(key)])
R_MM <- 1.35                                           # point radius in mm
split_disc <- function(fill, r = R_MM, n = 60) {
  k <- length(fill)
  do.call(grid::grobTree, lapply(seq_len(k), function(i) {
    th <- seq(pi/2 + (i-1)*2*pi/k, pi/2 + i*2*pi/k, length.out = n)
    grid::polygonGrob(x = unit(0.5, "npc") + unit(c(0, r*cos(th)), "mm"),
                      y = unit(0.5, "npc") + unit(c(0, r*sin(th)), "mm"),
                      gp = grid::gpar(fill = fill[i], col = NA)) }))
}

## plot
p <- ggplot() +
  geom_vline(xintercept=0.5, linetype=2, color="grey55") +
  geom_vline(xintercept=mP, linetype=3, color=CPRNP, linewidth=0.6) +
  geom_vline(xintercept=mG, linetype=3, color=CGEBV, linewidth=0.6) +
  geom_vline(xintercept=mR, linetype=3, color=CREM,  linewidth=0.6) +
  geom_segment(data=seg, aes(x=xmin, xend=xmax, y=y, yend=y), color="grey78", linewidth=1) +
  geom_point(data=filter(dl, !tied), aes(auc, y, color=pred), size=2.9) +
  scale_y_continuous(breaks=d$y, labels=d$lab, expand=expansion(add=0.6)) +
  scale_color_manual(values=cols) +
  scale_x_continuous(breaks=seq(0.4,0.9,0.1)) +
  guides(color=guide_legend(nrow=3)) +
  labs(x="out-of-farm prediction accuracy (AUC)", y=NULL) +
  TH + theme(legend.position="top", legend.justification="left", legend.text=element_markdown(size=12.5))

for (kk in unique(dl$key[dl$tied])) {
  s <- dl %>% filter(key == kk) %>% arrange(pred)
  p <- p + annotation_custom(split_disc(cols[as.character(s$pred)]),
                             xmin = s$auc[1], xmax = s$auc[1], ymin = s$y[1], ymax = s$y[1])
}

## save
ggsave("Figures/fig3_prnp_vs_gebv.png", p, width=9.5, height=6.2, dpi=170, bg="white", device=ragg::agg_png)
cat(sprintf("means: PRNP=%.2f GEBV=%.2f remainder=%.2f  (PRNP higher than GEBV in %d/13 farms)\nsaved fig3_prnp_vs_gebv.png\n",
    mP, mG, mR, sum(d$auc_PRNP > d$auc_GEBVfull)))
