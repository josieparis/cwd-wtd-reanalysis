# Step 1c: is the phenotype confounded with structure?
# binary (logistic) and tissue score (proportional odds) against region and the first 10 genomic PCs,
# i.e. how much of the phenotype structure alone explains

library(MASS)
library(dplyr)

## working dir
ROOT <- "<dir>/Seabury_data/reanalysis"

## read in PCs and phenotypes
d <- read.delim(file.path(ROOT,"01_structure","pca_merged.tsv")) %>%
  mutate(region = relevel(factor(region), ref="South"), sex = factor(sex_raw))
pcs <- paste0("PC",1:10)

## McFadden and Nagelkerke pseudo R2
mcf  <- function(m,m0) as.numeric(1 - logLik(m)/logLik(m0))
nagel<- function(m,m0){ n<-nobs(m); l0<-as.numeric(logLik(m0)); l1<-as.numeric(logLik(m))
  (1-exp((2/n)*(l0-l1)))/(1-exp((2/n)*l0)) }

## binary (logistic regression)
cat("=====================================================================\n")
cat("BINARY  CWD_Binary  (logistic regression)\n")
cat("=====================================================================\n")
b0  <- glm(binary~1, family=binomial, data=d)
bR  <- glm(binary~region, family=binomial, data=d)
bP  <- glm(reformulate(pcs,"binary"), family=binomial, data=d)
bP1 <- glm(binary~PC1, family=binomial, data=d)
bRP <- glm(reformulate(c("region",pcs),"binary"), family=binomial, data=d)
bCov<- glm(reformulate(c("region",pcs,"sex","age"),"binary"), family=binomial, data=d)
show <- function(tag,m){ cat(sprintf("  %-30s McFadden R2=%.3f  Nagelkerke=%.3f  AIC=%.1f\n",
        tag, mcf(m,b0), nagel(m,b0), AIC(m))) }
show("~ region (3 levels)", bR)
show("~ PC1 only", bP1)
show("~ PC1..PC10", bP)
show("~ region + PC1..PC10", bRP)
show("~ region + PCs + sex + age", bCov)
cat(sprintf("\n  LRT region vs null:      p=%.3e\n", anova(b0,bR,test="Chisq")$`Pr(>Chi)`[2]))
cat(sprintf("  LRT PC1..PC10 vs null:   p=%.3e\n", anova(b0,bP,test="Chisq")$`Pr(>Chi)`[2]))
cat(sprintf("  LRT PCs add over region: p=%.3e  (structure beyond region)\n",
    anova(bR,bRP,test="Chisq")$`Pr(>Chi)`[2]))

## tissue score 0/1/2 (proportional odds, polr)
cat("\n=====================================================================\n")
cat("ORDINAL  CMS_CWDScore 0/1/2  (proportional-odds logistic, polr)\n")
cat("=====================================================================\n")
d <- d %>% mutate(score = factor(score, ordered=TRUE))
o0 <- polr(score~1, data=d, Hess=TRUE)
oR <- polr(score~region, data=d, Hess=TRUE)
oP <- polr(reformulate(pcs,"score"), data=d, Hess=TRUE)
oRP<- polr(reformulate(c("region",pcs),"score"), data=d, Hess=TRUE)
showo<- function(tag,m){ cat(sprintf("  %-30s McFadden R2=%.3f  AIC=%.1f\n", tag, mcf(m,o0), AIC(m))) }
showo("~ region", oR)
showo("~ PC1..PC10", oP)
showo("~ region + PC1..PC10", oRP)

## region odds ratios (binary, South = reference)
cat("\n=====================================================================\n")
cat("Region coefficients (binary logistic, ref = South)  -- odds ratios\n")
cat("=====================================================================\n")
co <- summary(bR)$coefficients
for(r in rownames(co)) cat(sprintf("  %-18s OR=%6.2f  p=%.3e\n", r, exp(co[r,1]), co[r,4]))
cat("\nDONE\n")
