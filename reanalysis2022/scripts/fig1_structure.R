# Figure 1: CWD status is confounded with region, farm and relatedness
# a PCA of the GRM coloured by region, b GRM relatedness distribution, c CWD prevalence by farm (2022)
# GRM = PLINK 1.9 --make-rel (VanRaden) in both datasets; close relative = GRM >= 0.2, first degree = GRM >= 0.4

library(ggplot2)
library(patchwork)
library(dplyr)
library(ggrepel)

## working dir
setwd("<dir>/Seabury_data")

## colours, helpers and theme (region "East" is shown as Northeast)
regcol  <- c(Northeast="#8E6BB1", Midwest="#F2A03E", South="#1CA089")
fillcol <- setNames(scales::alpha(regcol, 0.8), names(regcol))
statcol <- c("CWD-"="grey60", "CWD+"="#B2182B")                 # control / case
relabR <- function(r) ifelse(r=="East","Northeast",r)
relabF <- function(f) gsub("East","Northeast",f)
sq     <- function(f) as.matrix(read.table(gzfile(f), sep="\t"))
crel   <- function(K,thr){ diag(K)<-0; 100*mean(rowSums(K>=thr)>=1) }
flipS  <- function(x,reg){ if(mean(x[reg=="South"])<0) -x else x }
update_geom_defaults("text", list(family="Arial")); update_geom_defaults(ggrepel::GeomTextRepel, list(family="Arial"))
th <- theme_classic(base_size=16, base_family="Arial") +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=16),
        plot.tag=element_text(face="bold", size=24), legend.position="none",
        axis.text=element_text(colour="black"))

## a: PCA and PC1 density by CWD status
p20 <- read.delim("reanalysis/01_structure/pca_merged.tsv", check.names=FALSE, stringsAsFactors=FALSE)
ev20 <- scan("reanalysis/01_structure/pca.eigenval", quiet=TRUE); v20 <- 100*ev20/sum(ev20[1:20])
p20$region <- factor(relabR(p20$region), levels=names(regcol))
p20$PC1 <- flipS(p20$PC1, as.character(p20$region))
p20$binary <- suppressWarnings(as.numeric(p20$binary))
pe <- read.table("reanalysis2022/01_structure/pca2022.eigenvec", stringsAsFactors=FALSE)
ev22 <- scan("reanalysis2022/01_structure/pca2022.eigenval", quiet=TRUE); v22 <- 100*ev22/sum(ev22[1:20])
cov <- read.csv("reanalysis2022/data/cov_1170.csv", stringsAsFactors=FALSE); cov$binary <- as.numeric(cov$binary)
r22 <- relabR(cov$region[match(pe$V2, cov$CMSID)])
d22 <- data.frame(PC1=pe$V3, PC2=pe$V4, region=factor(r22, levels=names(regcol)),
                  binary=as.numeric(cov$binary[match(pe$V2, cov$CMSID)]))
d22$PC1 <- flipS(d22$PC1, as.character(d22$region))

## same PC1 range for the scatter and the density plot, with room for the region labels
xrng <- function(df){ g<-mean(df$PC1); e<-g+(tapply(df$PC1,df$region,median)-g)*2.2
  r<-range(c(df$PC1,e),na.rm=TRUE); r+c(-1,1)*diff(r)*0.08 }
xr20 <- xrng(p20); xr22 <- xrng(d22)

pca_scatter <- function(df, v, ttl, xr){
  gc1 <- mean(df$PC1); gc2 <- mean(df$PC2)                       # push labels away from the centre
  cent <- df %>% group_by(region) %>% summarise(cx=median(PC1), cy=median(PC2), .groups="drop") %>%
    mutate(PC1=gc1+(cx-gc1)*2.2, PC2=gc2+(cy-gc2)*2.2)
  ggplot(df, aes(PC1, PC2)) +
    geom_point(aes(fill=region, color=region), shape=21, size=2, stroke=0.4) +
    geom_text_repel(data=cent, aes(label=region, color=region), fontface=2, size=6, seed=1,
                    box.padding=0.3, min.segment.length=Inf, bg.color="white", bg.r=0.2) +
    scale_fill_manual(values=fillcol) + scale_color_manual(values=regcol) +
    scale_y_continuous(expand=expansion(mult=0.12)) +
    coord_cartesian(xlim=xr, clip="off") +
    labs(x=NULL, y=sprintf("PC2 (%.1f%%)", v[2]), title=ttl) + th +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
}
pc1_density <- function(df, v, xr){                              # PC1 by CWD status
  df <- df[!is.na(df$binary),]
  df$status <- factor(ifelse(df$binary==1,"CWD+","CWD-"), levels=c("CWD-","CWD+"))
  rv  <- cor(df$PC1, df$binary)                                 # point-biserial r
  mus <- df %>% group_by(status) %>% summarise(m=mean(PC1), .groups="drop")
  ggplot(df, aes(PC1, fill=status, color=status)) +
    geom_density(alpha=0.45, linewidth=0.5) +
    geom_vline(data=mus, aes(xintercept=m, color=status), linetype=2, linewidth=0.5, show.legend=FALSE) +
    annotate("text", x=xr[2], y=Inf, hjust=1, vjust=1.5, size=4.4, fontface=3, color="grey20",
             label=sprintf("r = %.2f", rv)) +
    scale_fill_manual(values=statcol, name=NULL) + scale_color_manual(values=statcol, name=NULL) +
    scale_y_continuous(expand=expansion(mult=c(0,0.10))) +
    coord_cartesian(xlim=xr) +
    labs(x=sprintf("PC1 (%.1f%%)", v[1]), y="density") + th
}
pa20 <- pca_scatter(p20, v20, "2020 dataset (n = 807)",   xr20) + labs(tag="a")
pa22 <- pca_scatter(d22, v22, "2022 dataset (n = 1,170)", xr22)
pd20 <- pc1_density(p20, v20, xr20) +                           # CWD status legend on the left only
  theme(legend.position=c(0.015,0.96), legend.justification=c(0,1),
        legend.key.size=unit(0.85,"lines"), legend.text=element_text(size=11),
        legend.background=element_rect(fill="white", color=NA))
pd22 <- pc1_density(d22, v22, xr22)

## b: relatedness
G20 <- sq("reanalysis/01_structure/grm.rel.gz");            off20 <- G20[lower.tri(G20)]
G22 <- sq("reanalysis2022/01_structure/grm2022.rel.gz");    off22 <- G22[lower.tri(G22)]
fd20 <- sum(off20>=0.4); cr20 <- crel(G20,0.2); fd22 <- sum(off22>=0.4); cr22 <- crel(G22,0.2)
BRK <- seq(-0.20, 1.12, by=0.01)
mkh <- function(off){ h<-hist(off, breaks=BRK, plot=FALSE); data.frame(mid=h$mids, count=h$counts) }
h20 <- mkh(off20); h22 <- mkh(off22); YMAX <- max(h20$count, h22$count)*1.4
rel_panel <- function(hd, fd, cr, ttl){
  ggplot(hd[hd$count>0,], aes(mid, count)) +
    annotate("rect", xmin=0.4, xmax=1.10, ymin=0.5, ymax=Inf, fill="#fdd0d0", alpha=0.5) +
    geom_segment(aes(xend=mid, y=0.5, yend=count), color="grey35", linewidth=0.5) +
    geom_vline(xintercept=0.2, linetype=2, color="#e08214", linewidth=0.7) +
    geom_vline(xintercept=0.4, linetype=2, color="#b2182b", linewidth=0.7) +
    scale_y_log10(labels=scales::label_comma(), expand=expansion(mult=c(0,0))) +
    coord_cartesian(xlim=c(-0.12, 1.10), ylim=c(0.5, YMAX)) +
    labs(x="GRM off-diagonal (~2 x kinship)", y="pair count (log scale)", title=ttl) + th
}
pb20 <- rel_panel(h20, fd20, cr20, "2020 dataset relatedness") + labs(tag="b")
## 2022: colour pairs by within vs between farm (no farm IDs for 2020)
id22   <- read.table("reanalysis2022/01_structure/grm2022.rel.id", stringsAsFactors=FALSE)$V2
farm22 <- cov$farm[match(id22, cov$CMSID)]
same   <- outer(farm22, farm22, "==")[lower.tri(G22)]; same[is.na(same)] <- FALSE
hw <- mkh(off22[ same]); hb <- mkh(off22[!same])            # within / between counts per bin
h22c <- data.frame(mid=hw$mid, count=hw$count + hb$count,   # one bar per bin, coloured by the commoner type
                   cat=ifelse(hw$count > hb$count, "within farm", "between farm"))
h22c$cat <- factor(h22c$cat, levels=c("between farm","within farm"))
FARMCOL <- c("between farm"="grey70", "within farm"="#2166AC")
pb22 <- ggplot() +
  annotate("rect", xmin=0.4, xmax=1.10, ymin=0.5, ymax=Inf, fill="#fdd0d0", alpha=0.5) +
  geom_segment(data=subset(h22c, count>0),
               aes(x=mid, xend=mid, y=0.5, yend=count, color=cat), linewidth=0.6) +
  geom_vline(xintercept=0.2, linetype=2, color="#e08214", linewidth=0.7) +
  geom_vline(xintercept=0.4, linetype=2, color="#b2182b", linewidth=0.7) +
  scale_color_manual(values=FARMCOL, name=NULL) +
  scale_y_log10(labels=scales::label_comma(), expand=expansion(mult=c(0,0))) +
  coord_cartesian(xlim=c(-0.12, 1.10), ylim=c(0.5, YMAX)) +
  labs(x="GRM off-diagonal (~2 x kinship)", y="pair count (log scale)", title="2022 dataset relatedness") +
  th + theme(legend.position=c(0.955,0.999), legend.justification=c(1,1),  # lined up with the pink band edge
             legend.text=element_text(size=13), legend.key.height=unit(0.95,"lines"),
             legend.background=element_rect(fill="white", color="grey75", linewidth=0.3),
             legend.margin=margin(5,8,5,6))

## c: prevalence by farm
cov$region <- relabR(cov$region); cov$farm <- relabF(cov$farm)
cov2 <- cov %>% filter(!is.na(binary), !is.na(farm), farm!="")
rmean <- cov2 %>% group_by(region) %>% summarise(m=mean(binary)*100, .groups="drop")
lev <- rmean$region[order(rmean$m)]
ag <- cov2 %>% group_by(farm, region) %>% summarise(n=n(), prev=mean(binary)*100, .groups="drop") %>%
      filter(n >= 10) %>%                                          # the 16 farms with >= 10 deer
      mutate(region=factor(region, levels=lev)) %>% arrange(region, prev)
ag$farm <- factor(ag$farm, levels=ag$farm)
rmean$region <- factor(rmean$region, levels=lev)
lab <- setNames(sprintf("%s (region mean %.0f%%)", rmean$region, rmean$m), as.character(rmean$region))
pc <- ggplot(ag, aes(farm, prev, fill=region)) +
  geom_hline(data=rmean, aes(yintercept=m, color=region), linetype=2, linewidth=0.7, show.legend=FALSE) +
  geom_col(aes(color=region), linewidth=0.3, width=0.82) +
  geom_text(aes(label=paste0("n=", n)), angle=90, hjust=-0.15, size=4.3, color="grey25") +
  scale_fill_manual(values=fillcol, labels=lab, name=NULL) +
  scale_color_manual(values=regcol, guide="none") +
  scale_y_continuous(limits=c(0,123), expand=expansion(mult=c(0,0.02))) +
  labs(x=NULL, y="Farm CWD prevalence (%)", tag="c") +
  theme_classic(base_size=16, base_family="Arial") +
  theme(plot.tag=element_text(face="bold", size=24), axis.text=element_text(colour="black"),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=13),
        legend.text=element_text(size=14),
        legend.position="top", legend.justification="left", legend.direction="horizontal",
        legend.key.size=unit(1.2,"lines"), legend.margin=margin(0,0,0,0))

## assemble and save
fig <- (pa20 | pa22) / (pd20 | pd22) / (pb20 | pb22) / pc +
  plot_layout(heights=c(1, 0.34, 0.92, 1.15))
ggsave("Figures/fig1_structure.png", fig, width=12.5, height=16.5, dpi=170, bg="white", device=ragg::agg_png)
cat(sprintf("2020 PC1=%.1f%% fd=%d cr=%.0f%% | 2022 PC1=%.1f%% fd=%d cr=%.0f%% | saved fig1_structure.png\n",
            v20[1], fd20, cr20, v22[1], fd22, cr22))
