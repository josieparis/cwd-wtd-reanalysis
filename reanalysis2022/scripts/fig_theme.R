# Shared colours and ggplot2 theme for the figures, sourced by the figure scripts
# All text is Arial, so save with device = ragg::agg_png

library(ggplot2)
library(scales)

## colours
REGCOL  <- c(Northeast="#8E6BB1", Midwest="#F2A03E", South="#1CA089")          # regions
REGFILL <- setNames(scales::alpha(REGCOL, 0.8), names(REGCOL))
METHCOL <- c(naive="#b2182b", pc="#F2A03E", grm="#2166AC")                     # naive / +PCs / GRM
COHCOL  <- c(`2020`="#3B7DB3", `2022`="#C0392B")                               # 2020 / 2022 datasets
FONT    <- "Arial"

## theme
theme_seabury <- function(base=13)
  theme_classic(base_size=base, base_family=FONT) +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=base),
        plot.tag=element_text(face="bold", size=20),
        axis.text=element_text(size=rel(1), colour="black"),        # 13 pt
        axis.title=element_text(size=rel(15/13)),                   # 15 pt
        legend.text=element_text(size=rel(12.5/13)),                # 12.5 pt
        legend.title=element_blank(), legend.background=element_blank(),
        legend.key=element_blank())

## geoms do not inherit the theme font
update_geom_defaults("text",  list(family=FONT))
update_geom_defaults("label", list(family=FONT))
if (requireNamespace("ggrepel", quietly=TRUE)) {
  update_geom_defaults(ggrepel::GeomTextRepel,  list(family=FONT))
  update_geom_defaults(ggrepel::GeomLabelRepel, list(family=FONT))
}
