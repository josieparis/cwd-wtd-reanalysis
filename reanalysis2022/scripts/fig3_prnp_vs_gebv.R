# Figure 3: out-of-farm AUC per farm for PRNP codon 96 alone, the full GEBV
# and the GEBV without the PRNP chromosome (polygenic remainder)

library(ggplot2)
library(dplyr)
library(ggtext)

## working dir
setwd("<dir>/Seabury_data")

source("reanalysis2022/scripts/fig_theme.R"); TH <- theme_seabury()
CPRNP <- "#1B7837"; CGEBV <- "#B2182B"; CREM <- "#E08214"

## read in per-farm AUCs
d <- read.csv("reanalysis2022/data/prnp_threeway_perfarm.csv", stringsAsFactors=FALSE)
d$farm <- gsub("East","Northeast", d$farm)
d <- d[order(d$n),]                                              # order farms by sample size
d$star <- pmax(d$auc_GEBVfull, d$auc_GEBVminus13) > d$auc_PRNP   # asterisk if a GEBV beats PRNP
d$lab  <- sprintf("%s (n=%d)%s", d$farm, d$n, ifelse(d$star, "*", ""))
d$y    <- seq_len(nrow(d))
mP <- mean(d$auc_PRNP); mG <- mean(d$auc_GEBVfull); mR <- mean(d$auc_GEBVminus13)
lv <- c(sprintf("*PRNP* codon-96 alone (mean %.2f)", mP),
        sprintf("full genome-wide GEBV (mean %.2f)", mG),
        sprintf("GEBV minus *PRNP* chromosome (mean %.2f)", mR))
dl <- rbind(data.frame(y=d$y, auc=d$auc_PRNP,        pred=lv[1]),
            data.frame(y=d$y, auc=d$auc_GEBVfull,    pred=lv[2]),
            data.frame(y=d$y, auc=d$auc_GEBVminus13, pred=lv[3]))
dl$pred <- factor(dl$pred, levels=lv)
seg <- data.frame(y=d$y, xmin=pmin(d$auc_PRNP,d$auc_GEBVfull,d$auc_GEBVminus13),
                  xmax=pmax(d$auc_PRNP,d$auc_GEBVfull,d$auc_GEBVminus13))
cols <- setNames(c(CPRNP,CGEBV,CREM), lv)

## tied AUCs are drawn as one point split between the colours
dl$key <- paste(dl$y, dl$auc); tied <- dl$key %in% dl$key[duplicated(dl$key)]
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
  geom_point(data=dl[!tied, ], aes(auc, y, color=pred), size=2.9) +
  scale_y_continuous(breaks=d$y, labels=d$lab, expand=expansion(add=0.6)) +
  scale_color_manual(values=cols) +
  scale_x_continuous(breaks=seq(0.4,0.9,0.1)) +
  guides(color=guide_legend(nrow=3)) +
  labs(x="out-of-farm prediction accuracy (AUC)", y=NULL) +
  TH + theme(legend.position="top", legend.justification="left", legend.text=element_markdown(size=12.5))

for (kk in unique(dl$key[tied])) {
  s <- dl[dl$key == kk, ]; s <- s[order(s$pred), ]
  p <- p + annotation_custom(split_disc(cols[as.character(s$pred)]),
                             xmin = s$auc[1], xmax = s$auc[1], ymin = s$y[1], ymax = s$y[1])
}

## save
ggsave("Figures/fig3_prnp_vs_gebv.png", p, width=9.5, height=6.2, dpi=170, bg="white", device=ragg::agg_png)
cat(sprintf("means: PRNP=%.2f GEBV=%.2f remainder=%.2f  (PRNP higher than GEBV in %d/13 farms)\nsaved fig3_prnp_vs_gebv.png\n",
    mP, mG, mR, sum(d$auc_PRNP > d$auc_GEBVfull)))
