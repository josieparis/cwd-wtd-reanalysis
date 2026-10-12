# Supplementary Figure 2: GWAS calibration, only PRNP survives structure correction (binary CWD, both datasets)
# a QQ plots (naive, +10 PCs, single GRM), b Manhattan plots for the single-GRM model (GCTA --mlma, as in Table 1),
# c number of SNPs at P < 5e-5 per model, d tissue score as a covariate (2020, GEMMA)

library(tidyverse)
library(patchwork)

## working dir
setwd("<dir>/Seabury_data")

source("reanalysis2022/scripts/fig_theme.R")
PRNP <- "G96S_G_to_A_at286"; TH <- theme_seabury()
## genomic control lambda
lam <- function(p){ p<-p[is.finite(p)&p>0]; median(qchisq(p,1,lower.tail=FALSE))/qchisq(0.5,1) }

## read in data
bim <- read.table("reanalysis/00_data/deer.bim", stringsAsFactors=FALSE,
                  col.names=c("chr","marker","cm","pos","a1","a2"))
m20 <- read.delim("reanalysis/02_assoc/output/gcta_mlma_bin.mlma", stringsAsFactors=FALSE)   # GCTA --mlma (Table 1)
g20 <- read.delim(gzfile("reanalysis/01_structure/gwas_binary.tsv.gz"), stringsAsFactors=FALSE) %>%
  inner_join(bim %>% select(marker, chr, pos), by="marker") %>% arrange(marker) %>%
  left_join(m20 %>% select(marker=SNP, p_grm=p), by="marker")
stopifnot(!anyNA(g20$p_grm))
g22 <- read.delim(gzfile("reanalysis2022/data/gwas2022_binary.tsv.gz"), stringsAsFactors=FALSE)
gem <- function(f) read.delim(file.path("reanalysis/02_assoc/output", f), stringsAsFactors=FALSE) %>% select(marker=rs, p_wald)
cl  <- inner_join(gem("gwaa_bin_cov.assoc.txt")      %>% rename(p_base=p_wald),     # binary ~ sex/age/region
                  gem("gwaa_bin_covscore.assoc.txt") %>% rename(p_score=p_wald),    # + tissue score
                  by="marker") %>% arrange(marker)
auto <- read.table("reanalysis/02_assoc/deer_bin.bim", stringsAsFactors=FALSE)$V2  # gwaa_* scans include chrX
cl  <- cl %>% filter(marker %in% auto)                                             # keep the Table 1 autosomal set

## a: QQ plots
qq_thin <- function(p){
  p <- sort(p[is.finite(p)&p>0]); n <- length(p)
  data.frame(exp=-log10(ppoints(n)), obs=-log10(p)) %>%
    filter(obs>1.5 | row_number() %% 40 == 0)                  # keep tail, thin the null bulk
}
modname <- c(naive="naive", pc="+10 PCs", grm="single GRM")
qq_panel <- function(pl, ttl){
  labs <- vapply(names(pl), function(m) sprintf("%s (λ=%.2f)", modname[m], lam(pl[[m]])), "")
  d <- map2_dfr(names(pl), labs, function(m,l) qq_thin(pl[[m]]) %>% mutate(model=l)) %>%
    mutate(model=factor(model, levels=labs))
  ggplot(d, aes(exp, obs, color=model)) +
    geom_abline(slope=1, intercept=0, linetype=2, color="grey55") +
    geom_hline(yintercept=-log10(5e-8), linetype=3, color="grey70") +
    geom_point(size=0.5) +
    scale_color_manual(values=setNames(unname(METHCOL), labs)) +
    guides(color=guide_legend(override.aes=list(size=2.5))) +
    labs(x=expression(expected~-log[10](italic(p))), y=expression(observed~-log[10](italic(p))), title=ttl) +
    TH + theme(legend.position=c(0.03,0.97), legend.justification=c(0,1), legend.text=element_text(size=11.5))
}
qa <- qq_panel(list(naive=g20$p_naive, pc=g20$p_pc, grm=g20$p_grm), "2020 dataset") + labs(tag="a")
qb <- qq_panel(list(naive=g22$p_naive, pc=g22$p_pc, grm=g22$p_grm),   "2022 dataset")

## b: Manhattan plots (single GRM)
manh <- function(df, pcol, ttl){
  d <- df %>% filter(is.finite(.data[[pcol]]), !is.na(chr)) %>%
    mutate(chr=as.integer(chr), bp0=if("pos" %in% names(df)) pos else bp) %>%
    arrange(chr, bp0)
  ax <- d %>% group_by(chr) %>% summarise(len=max(bp0)) %>%       # chromosome offsets and label positions
    mutate(off=cumsum(c(0, head(len, -1))), center=off + len/2)
  d <- d %>% left_join(ax %>% select(chr, off), by="chr") %>% mutate(x=bp0 + off)
  prnp <- d %>% filter(marker==PRNP)
  ggplot(d, aes(x, -log10(.data[[pcol]]))) +
    geom_hline(yintercept=-log10(5e-8), linetype=2, color="#b2182b", linewidth=0.4) +
    geom_hline(yintercept=-log10(5e-5), linetype=3, color="#F2A03E", linewidth=0.4) +
    geom_point(aes(color=factor(chr %% 2)), size=0.45) +
    scale_color_manual(values=c("0"="grey72","1"="grey45"), guide="none") +
    geom_point(data=prnp, color="#b2182b", size=2.6, shape=18) +
    geom_text(data=prnp, aes(label="PRNP"), vjust=-0.9, color="#b2182b", fontface=3, size=3.4) +
    scale_x_continuous(breaks=ax$center, labels=ax$chr, expand=expansion(mult=0.01)) +
    scale_y_continuous(expand=expansion(mult=c(0,0.12)), limits=c(0, NA)) +
    labs(x="chromosome", y=expression(-log[10](italic(p))), title=ttl) +
    TH + theme(panel.grid=element_blank(), axis.text.x=element_text(size=9.5))
}
mb20 <- manh(g20, "p_grm",   "2020 dataset, single-GRM model") + labs(tag="b")
mb22 <- manh(g22, "p_grm",   "2022 dataset, single-GRM model")

## c: candidate SNP counts (P < 5e-5)
cnt <- function(p) sum(p < 5e-5, na.rm=TRUE)
cc <- data.frame(
  cohort = rep(c("2020","2022"), each=3),
  model  = factor(rep(c("naive","+10 PCs","single GRM"),2), levels=c("naive","+10 PCs","single GRM")),
  count  = c(cnt(g20$p_naive), cnt(g20$p_pc), cnt(g20$p_grm),
             cnt(g22$p_naive), cnt(g22$p_pc), cnt(g22$p_grm)))
cc <- cc %>% group_by(model) %>% mutate(tie = n_distinct(count) == 1, hi = count == max(count)) %>% ungroup() %>%
  mutate(hj = if_else(model == "naive", 1.2, 0.5),                   # naive labels to the left of the lines
         vj = case_when(model == "naive" & hi ~ 0.1, model == "naive" ~ 0.9, hi ~ -0.9, TRUE ~ 1.9))
pc <- ggplot(cc, aes(model, count, color=cohort, group=cohort)) +
  geom_line(linewidth=0.8) + geom_point(size=3) +
  geom_text(data=filter(cc, !tie), aes(label=scales::comma(count), hjust=hj, vjust=vj), size=4, show.legend=FALSE) +
  geom_text(data=filter(cc, tie, cohort=="2020"), aes(label=scales::comma(count)), vjust=-0.9, size=4, color="black") +
  scale_color_manual(values=COHCOL) +
  scale_y_log10(labels=scales::label_comma(), expand=expansion(mult=c(0.08,0.18))) +
  labs(x=NULL, y=expression(SNPs~italic(P)<5%*%10^-5~(log)), tag="c") +
  TH + theme(legend.position=c(0.97,0.97), legend.justification=c(1,1))

## d: tissue score as a covariate (2020, GEMMA scans with and without the score)
## the score contains the outcome (binary = score > 0), so PRNP is lost and a chr21 SNP goes up
cl2 <- cl %>% filter(is.finite(p_base), is.finite(p_score))
prnp <- cl2 %>% filter(marker==PRNP); chr21 <- cl2 %>% arrange(p_score) %>% slice(1)
stopifnot(chr21$marker == "Affx-573993815")                      # the chr21 SNP named in Results
pd <- ggplot(cl2, aes(-log10(p_base), -log10(p_score))) +
  geom_abline(slope=1, linetype=2, color="grey55") +
  geom_hline(yintercept=-log10(5e-8), linetype=3, color="grey70") +
  geom_vline(xintercept=-log10(5e-8), linetype=3, color="grey70") +
  geom_point(size=0.5, color="grey72", alpha=0.5) +
  geom_point(data=prnp,  color="#b2182b", size=3) +
  geom_point(data=chr21, color="#2166AC", size=3) +
  geom_text(data=prnp,  aes(label="PRNP\nlost with score"), color="#b2182b", size=3.9, hjust=1.1, vjust=0.3) +
  geom_text(data=chr21, aes(label="chr21 SNP, promoted by score"), color="#2166AC", size=3.7, hjust=-0.04, vjust=0.5) +
  annotate("text", x=Inf, y=4.4, hjust=1.05, size=3.7, color="grey35",
           label=sprintf("λ = %.2f without score\nλ = %.2f with score", lam(cl2$p_base), lam(cl2$p_score))) +
  labs(x=expression(-log[10](italic(p))*",  binary ~ sex + age + region"),
       y=expression(-log[10](italic(p))*",  ... + tissue score"), tag="d") + TH

## assemble and save
fig <- (qa | qb) / mb20 / mb22 / (pc | pd) + plot_layout(heights=c(1.05, 0.62, 0.62, 1.0))
ggsave("Figures/figS2_gwas_calibration.png", fig, width=12.5, height=15.5, dpi=170, bg="white", device=ragg::agg_png)
cat(sprintf("QQ lambda 2020 naive=%.1f pc=%.2f grm=%.2f | 2022 naive=%.1f pc=%.2f grm=%.2f\n",
    lam(g20$p_naive),lam(g20$p_pc),lam(g20$p_grm), lam(g22$p_naive),lam(g22$p_pc),lam(g22$p_grm)))
cat(sprintf("counts: %s\n", paste(cc$cohort, cc$model, cc$count, collapse=" | ")))
cat(sprintf("panel d (GEMMA): PRNP %.1e -> %.3f | %s %.1e -> %.1e | lambda %.2f -> %.2f\n",
    prnp$p_base, prnp$p_score, chr21$marker, chr21$p_base, chr21$p_score, lam(cl2$p_base), lam(cl2$p_score)))
cat("saved figS2_gwas_calibration.png\n")
