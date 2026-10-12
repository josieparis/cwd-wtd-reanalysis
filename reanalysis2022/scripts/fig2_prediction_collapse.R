# Figure 2: genomic prediction does not transfer to unseen farms (2022 dataset)
# a per-farm AUC, size-matched CV vs leave-one-farm-out, b AUC by relatedness to the training set,
# c all metrics for the three validation designs, d within vs across region prediction

library(ggplot2)
library(patchwork)
library(dplyr)

## working dir
setwd("<dir>/Seabury_data")

source("reanalysis2022/scripts/fig_theme.R"); TH <- theme_seabury()        # Arial, 13/15 pt text
## mean and SD from "mean (SD)" strings
num <- function(x) as.numeric(sub(" .*","",x))
sdv <- function(x) as.numeric(gsub(".*\\((.*)\\).*","\\1", x))
mlab <- c(auc="AUC", acc="Classification\naccuracy", sens="Sensitivity", spec="Specificity", mcc="MCC", r="Pred. r", rmse="RMSE")
COLL <- "#b2182b"; CMATCH <- "#3B7DB3"; CFOLD <- "grey55"
## prediction_lofo.R results (PLINK GRM, no covariates, the Table 2 run), used in a and c
PR <- readRDS("reanalysis2022/data/prediction_lofo_plink.rds")$none

## a: per-farm AUC
cv <- read.csv("reanalysis2022/data/cov_1170.csv", stringsAsFactors=FALSE)
ycv <- suppressWarnings(as.numeric(cv$binary))
pf <- data.frame(farm=rownames(PR$lofo), auc_lofo=PR$lofo[, "auc"], match=PR$size_matched[, "auc"], stringsAsFactors=FALSE)
pf$n <- sapply(pf$farm, function(f) sum(cv$farm == f & !is.na(ycv)))
pf$farm <- gsub("East","Northeast",pf$farm)
pf <- pf[order(pf$n),]; pf$lab <- factor(sprintf("%s (n=%d)", pf$farm, pf$n), levels=sprintf("%s (n=%d)", pf$farm, pf$n))
pfl <- rbind(data.frame(lab=pf$lab, auc=pf$auc_lofo, design="leave-one-farm-out"),
             data.frame(lab=pf$lab, auc=pf$match,    design="sample size-matched CV"))
pa <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=0.5, ymin=-Inf, ymax=Inf, fill="grey92") +
  geom_segment(data=pf, aes(x=auc_lofo, xend=match, y=lab, yend=lab), color="grey70", linewidth=1) +
  geom_vline(xintercept=0.5, linetype=2, color="grey45") +
  geom_point(data=pfl, aes(auc, lab, color=design), size=2.8) +
  scale_color_manual(values=c(`leave-one-farm-out`=COLL, `sample size-matched CV`=CMATCH)) +
  scale_x_continuous(limits=c(0.33,0.95), breaks=seq(0.4,0.9,0.1)) +
  labs(x="prediction accuracy (AUC)", y=NULL, tag="a") +
  TH + theme(legend.position="top", legend.justification="left")

## b: AUC by the test deer's closest relative in the training set (random 5-fold, relatedness_binned_auc.R)
rb <- read.csv("reanalysis2022/data/relatedness_binned_auc.csv", stringsAsFactors=FALSE)
rb$pct <- rb$n / sum(rb$n) * 100; rb$ndeer <- round(rb$n / 10)     # deer per CV replicate (10 replicates)
deg <- gsub("-", "–", sub(" deg$", " degree", rb$degree))           # e.g. 1st-2nd degree, with an en dash
xl <- sprintf("%s\n(%s)", deg, gsub("-", "–", rb$bin)); rb$xlab <- factor(xl, levels=xl)        # GRM range under each label
pb <- ggplot(rb, aes(xlab, auc)) +
  geom_hline(yintercept=0.5, linetype=2, color="grey45") +
  geom_col(fill=CMATCH, color=CMATCH, alpha=0.75, linewidth=0.5, width=0.7, na.rm=TRUE) +   # same style as c and d
  geom_text(aes(label=ifelse(is.na(auc),"",sprintf("AUC %.2f",auc))), vjust=-0.6, size=4.4, na.rm=TRUE) +
  geom_text(aes(y=0.03, label=sprintf("%s deer\n(%.0f%%)", format(ndeer, big.mark=","), pct)),
            vjust=0, size=4.1, lineheight=0.9, color="black") +
  scale_y_continuous(limits=c(0,0.93), expand=expansion(mult=c(0,0.04))) +
  labs(x="closest relative in the training set (GRM)",
       y="prediction accuracy (AUC)", tag="b") +
  TH + theme(axis.text.x=element_text(size=11, lineheight=0.95))

## c: all metrics for the three validation designs
des <- c(random_5fold="random 5-fold", size_matched="sample size-matched", lofo="leave-one-farm-out")
lc <- do.call(rbind, lapply(names(des), function(d) { m <- PR[[d]]          # folds (5-fold) or held-out farms (others)
  data.frame(metric=colnames(m), design=des[d], value=colMeans(m, na.rm=TRUE), sd=apply(m, 2, sd, na.rm=TRUE)) }))
lc$metric <- factor(mlab[lc$metric], levels=mlab); lc$design <- factor(lc$design, levels=des)
pc <- ggplot(lc, aes(metric, value, fill=design, color=design)) +
  geom_hline(yintercept=c(0,0.5), linetype=c(1,2), color=c("grey30","grey55")) +
  geom_col(position=position_dodge(0.8), width=0.72, linewidth=0.5, alpha=0.75) +   # border = bar colour, translucent fill
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), position=position_dodge(0.8), width=0.25, linewidth=0.3, color="black") +
  scale_fill_manual(values=setNames(c(CFOLD,CMATCH,COLL), des)) +
  scale_color_manual(values=setNames(c(CFOLD,CMATCH,COLL), des)) +
  labs(x=NULL, y="value (mean ± SD)", tag="c") +
  TH + theme(legend.position="top", legend.justification="left")

## d: within vs across region (within_vs_across_region.R)
wr <- read.csv("reanalysis2022/data/within_vs_across_region_summary.csv", stringsAsFactors=FALSE)
ser <- c(within_region="within region", across_region="across region", across_region_plusPC="across region + PCs")
ld <- do.call(rbind, lapply(names(ser), function(s)
  data.frame(metric=wr$metric, series=ser[s], value=num(wr[[s]]), sd=sdv(wr[[s]]))))
ld$metric <- factor(mlab[ld$metric], levels=mlab); ld$series <- factor(ld$series, levels=ser)
pd <- ggplot(ld, aes(metric, value, fill=series, color=series)) +
  geom_hline(yintercept=c(0,0.5), linetype=c(1,2), color=c("grey30","grey55")) +
  geom_col(position=position_dodge(0.8), width=0.72, linewidth=0.5, alpha=0.75) +   # border = bar colour, translucent fill
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), position=position_dodge(0.8), width=0.25, linewidth=0.3, color="black") +
  scale_fill_manual(values=setNames(c("#1CA089","#F2A03E",COLL), ser)) +
  scale_color_manual(values=setNames(c("#1CA089","#F2A03E",COLL), ser)) +
  labs(x=NULL, y="value (mean ± SD)", tag="d") +
  TH + theme(legend.position="top", legend.justification="left")

## assemble and save
fig <- ((pa | pb) + plot_layout(widths=c(0.95, 1.15))) / pc / pd + plot_layout(heights=c(1.05, 0.85, 0.85))
ggsave("Figures/fig2_prediction_collapse.png", fig, width=12.5, height=14.5, dpi=170, bg="white",
       device=ragg::agg_png)
cat("saved fig2_prediction_collapse.png\n")
