# Supplementary Figure 3: ARSB intron-4 minor allele frequency vs farm CWD prevalence in each region
# each regional correlation rests on a few farms and one dominant herd (ringed), so the gene x region
# interaction reflects which herds were sampled rather than a biological gradient

library(ggplot2)
library(patchwork)
library(ggrepel)
library(ggtext)
library(dplyr)

## working dir
setwd("<dir>/Seabury_data")

source("reanalysis2022/scripts/fig_theme.R"); TH <- theme_seabury()

## read in per-farm data (region "East" shown as Northeast)
ag <- read.csv("reanalysis2022/data/arsb_perfarm.csv", stringsAsFactors=FALSE) %>%
  mutate(region=gsub("East","Northeast", region), farm=gsub("East","Northeast", farm), prevpct=prev*100)
## dominant herd in each region, dropped for the second correlation
dropfarm <- c(Northeast="Northeast2", Midwest="Midwest6", South="South5")

## one panel per region
panel <- function(R, tag){
  s <- ag %>% filter(region==R); col <- REGCOL[[R]]
  ct <- cor.test(s$prev, s$maf); r <- unname(ct$estimate); p <- ct$p.value
  df <- dropfarm[[R]]; keep <- s %>% filter(farm!=df); rdrop <- cor(keep$prev, keep$maf)
  ggplot(s, aes(prevpct, maf)) +
    geom_smooth(method="lm", formula=y~x, se=FALSE, color=col, linetype=2, linewidth=0.7) +
    geom_point(aes(size=n), shape=21, fill=scales::alpha(col,0.8), color="grey30", stroke=0.4) +
    geom_point(data=filter(s, farm==df), shape=1, size=7, color="grey15", stroke=1.1) +
    geom_text_repel(aes(label=farm), size=3.8, seed=1, box.padding=0.6, point.padding=0.4,
                    min.segment.length=0, segment.color="grey65", segment.size=0.3,
                    bg.color="white", bg.r=0.15, max.overlaps=Inf) +
    scale_size_area(max_size=11, limits=c(0, max(ag$n)), breaks=c(25,100,300), name="farm n") +
    scale_x_continuous(expand=expansion(mult=0.13)) + scale_y_continuous(expand=expansion(mult=0.12)) +
    guides(size=guide_legend(override.aes=list(fill="grey60", color="grey35"))) +
    labs(x="farm CWD prevalence (%)", y="*ARSB* intron-4 minor-allele frequency",
         title=sprintf("%s (%d farms)", R, nrow(s)),
         subtitle=sprintf("*r* = %+.2f, *p* = %.3f &nbsp;&nbsp;(drop %s: *r* = %+.2f)", r, p, df, rdrop),
         tag=tag) +
    TH + theme(axis.title.y=element_markdown(), plot.subtitle=element_markdown(size=11.5, hjust=0.5))
}

## plot and save
fig <- (panel("Northeast","a") | panel("Midwest","b") | panel("South","c")) +
  plot_layout(guides="collect") & theme(legend.position="bottom")
ggsave("Figures/figS3_arsb.png", fig, width=15, height=5.8, dpi=170, bg="white", device=ragg::agg_png)
cat("saved Figures/figS3_arsb.png\n")
