# Figure 4: CWD heritability in the 2022 dataset is inflated by structure and not identifiable
# a observed-scale h2 for each covariate set and tool (Table 3, from build_table3.R, plus GxEMM/LDAK)
# b liability-scale h2 against the assumed population prevalence K (Lee transform of the GCTA estimates)

library(ggplot2)
library(patchwork)
library(dplyr)

## working dir
setwd("<dir>/Seabury_data")

source("reanalysis2022/scripts/fig_theme.R"); TH <- theme_seabury()
TOOLCOL <- c(`original value (Seabury et al. 2022)`="#111111", GCTA="#1F6F3F", GEMMA="#2C7FB8",
             rrBLUP="#E08A1E", `GxEMM/LDAK`="#8856A7")

## a: observed-scale h2 by covariate set and tool
cs <- c("no covariates","+ sex/age/region","+ sex/age/farm")
tl <- names(TOOLCOL)
t3 <- read.csv("MANUSCRIPT_tables/table3_heritability_observed_long.csv", stringsAsFactors=FALSE)  # from build_table3.R
dA <- t3 %>%
  filter(cohort == 2022) %>%
  transmute(covset = sub("^none$", "no covariates", covset),
            tool   = sub("Original (EMMAX)", "original value (Seabury et al. 2022)", tool, fixed=TRUE),
            h2, se)
## GxEMM/LDAK: HOM model, observed scale, 1,120 deer on farms with >= 10 deer and both cases and controls
fa <- readRDS("reanalysis2022/cluster_gxemm/gxemm_farmadj_results.rds")
gx <- function(m) { o <- fa[[m]]$obs$hom; c(unname(o$h2[1]), sqrt(o$h2Covmat[1, 1])) }
gR <- gx("A_region_adj"); gF <- gx("B_farm_adj")                # GxEMM/LDAK, region- and farm-adjusted
dA <- dA %>% bind_rows(data.frame(covset=cs[2:3], tool="GxEMM/LDAK", h2=c(gR[1], gF[1]), se=c(gR[2], gF[2])))
stopifnot(all(dA$covset %in% cs), all(dA$tool %in% tl))
dA <- dA %>% mutate(covset=factor(covset, levels=cs), tool=factor(tool, levels=tl))
pa <- ggplot(dA, aes(covset, h2, color=tool)) +
  geom_errorbar(aes(ymin=h2-se, ymax=h2+se), position=position_dodge(0.6), width=0.25, linewidth=0.5, na.rm=TRUE) +
  geom_point(position=position_dodge(0.6), size=2.6) +
  annotate("text", x=0.45, y=0.60, label="~0.60", color="grey35", size=4.2, fontface=2) +
  annotate("text", x=3.55, y=0.33, label="~0.32", color="grey35", size=4.2, fontface=2) +
  scale_color_manual(values=TOOLCOL) +
  scale_x_discrete(expand=expansion(add=0.75)) +
  scale_y_continuous(limits=c(0.20,0.73), breaks=seq(0.2,0.7,0.1)) +
  labs(x=NULL, y="observed-scale heritability (± SE)", tag="a") +
  TH + theme(legend.position=c(0.02,0.02), legend.justification=c(0,0), legend.text=element_text(size=11.5))

## b: liability-scale h2 vs prevalence
## no-covariate GCTA estimate (PLINK GRM); the points are the GCTA --prevalence values in Table 4
read_hsq <- function(f) { L <- strsplit(readLines(f), "\t")
  get <- function(key) as.numeric(Filter(function(z) z[1] == key, L)[[1]][-1])[1]
  list(obs = get("V(G)/Vp"), liab = get("V(G)/Vp_L")) }
KS  <- c(0.05, 0.10, 0.26)
hs  <- lapply(sprintf("reanalysis2022/cluster_h2/gcta_plink/results/h2_nocov_liab_%.2f.hsq", KS), read_hsq)
cov <- read.csv("reanalysis2022/data/cov_1170.csv", stringsAsFactors=FALSE)
P   <- mean(suppressWarnings(as.numeric(cov$binary))==1, na.rm=TRUE)   # sample case fraction 0.352
h2o <- hs[[1]]$obs                                                    # GCTA observed h2, no covariates (0.593)
liabf <- function(K){ Th<-qnorm(1-K); z<-dnorm(Th); h2o*(K*(1-K))^2/(P*(1-P)*z^2) }
dB <- data.frame(K=seq(0.01,0.50,0.002)) %>% mutate(h2L=liabf(K))
anchor <- data.frame(K=KS, h2L=sapply(hs, `[[`, "liab")) %>% mutate(lab=sprintf("%.2f", h2L))
lit    <- data.frame(K=c(0.34,0.42,0.49)) %>%
  mutate(h2L=liabf(K), lab=sprintf("%.2f", h2L),
         vj=if_else(h2L < 1, 2.2, -1.2))                                  # label below points under the h2 = 1 line
Kx <- uniroot(function(k) liabf(k) - 1, c(0.2, 0.5))$root                # where the curve crosses 1
## same transform for the covariate models (each with its own observed h2 and sample case fraction)
hs2 <- function(f){ L <- strsplit(readLines(f), "\t"); as.numeric(Filter(function(z) z[1]=="V(G)/Vp", L)[[1]][2]) }
c2  <- cov %>% mutate(cwd=suppressWarnings(as.numeric(binary)), sexn=suppressWarnings(as.numeric(sex)))
big2 <- c2 %>% filter(!is.na(cwd)) %>% count(farm) %>% filter(n >= 10) %>% pull(farm)
s1  <- c2 %>% filter(!is.na(cwd), !is.na(sexn), !is.na(age)); s2 <- s1 %>% filter(farm %in% big2)   # 1,166 / 1,151 deer
lee <- function(h, Pp, K){ z <- dnorm(qnorm(1-K)); h*(K*(1-K))^2/(Pp*(1-Pp)*z^2) }
GP  <- "reanalysis2022/cluster_h2/gcta_plink/results"
CS  <- c("no covariates", "+ sex/age/region", "+ sex/age/farm")
dC  <- bind_rows(data.frame(K=dB$K, h2L=dB$h2L, model=CS[1]),
                 data.frame(K=dB$K, h2L=lee(hs2(file.path(GP,"h2_region.hsq")), mean(s1$cwd), dB$K), model=CS[2]),
                 data.frame(K=dB$K, h2L=lee(hs2(file.path(GP,"h2_farm.hsq")),   mean(s2$cwd), dB$K), model=CS[3])) %>%
  mutate(model=factor(model, levels=CS))
pb <- ggplot(dB, aes(K, h2L)) +
  annotate("rect", xmin=0, xmax=0.5, ymin=1, ymax=Inf, fill="#F6D5D5") +       # Inf = panel edge
  geom_hline(yintercept=1, linetype=2, color="grey40") +
  geom_line(data=dC, aes(K, h2L, color=model, linetype=model, linewidth=model)) +
  scale_color_manual(values=setNames(c("#C0392B","grey35","grey35"), CS)) +
  scale_linetype_manual(values=setNames(c("solid","22","4212"), CS)) +
  scale_linewidth_manual(values=setNames(c(1.1,0.8,0.8), CS)) +
  geom_segment(x=Kx, xend=Kx, y=0.15, yend=1, linetype=3, color="#C0392B") +
  geom_point(data=anchor, color="#111111", size=2.4) +
  geom_text(data=anchor, aes(label=lab), color="#111111", size=4, hjust=1.35, vjust=-0.5) +
  geom_point(data=lit, shape=17, color="#24507A", size=2.6) +
  geom_text(data=lit, aes(label=lab, vjust=vj), color="#24507A", size=4) +
  annotate("text", x=0.32, y=1.07, label=sprintf("h² = 1 at K = %.2f", Kx), color="black", size=4.4, hjust=1) +
  scale_x_continuous(limits=c(0,0.5), expand=expansion(mult=c(0,0.01))) +
  scale_y_continuous(limits=c(0.15,1.15), breaks=seq(0.2,1.1,0.1), expand=expansion(mult=0)) +
  labs(x="assumed CWD population prevalence (K)", y="liability-scale heritability", tag="b") + TH +
  theme(legend.position=c(0.2,0.01), legend.justification=c(0,0), legend.key.width=unit(1.6,"lines"))   # under the farm curve

## assemble and save
fig <- pa | pb
ggsave("Figures/fig4_heritability.png", fig, width=12.5, height=5.6, dpi=170, bg="white", device=ragg::agg_png)
cat(sprintf("GxEMM/LDAK observed: region %.3f (SE %.3f) | farm %.3f (SE %.3f)\n", gR[1], gR[2], gF[1], gF[2]))
cat(sprintf("P(sample)=%.3f h2obs=%.3f | anchors %s | K=.34->%.2f .49->%.2f | crosses 1 at K=%.3f\n",
    P, h2o, paste(anchor$lab, collapse="/"), liabf(0.34), liabf(0.49), Kx))
cat("saved Figures/fig4_heritability.png\n")
