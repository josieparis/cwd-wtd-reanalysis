# Table 5: GxEMM GxE tests, and Supplementary Table 5: per-farm heritability, from the PLINK GRM GxEMM runs
# liability scale at K = 0.05, 0.10 and 0.26; region tests on all 1,166 complete-case deer, on the 1,120 deer
# of the farm models (Model A) and with farm as a fixed effect (Model B)

library(GxEMM)

## working dir
setwd("<dir>/Seabury_data")

OUT <- "MANUSCRIPT_tables"

## formatting helpers
se   <- function(C) sqrt(pmax(diag(as.matrix(C)), 0))
mns  <- function(s) gsub("-", "−", s)                      # hyphen to minus sign
num  <- function(x) if(!is.finite(x)) "NA" else if(abs(x) >= 100) sprintf("%.0f", x) else sprintf("%.2f", x)
pm   <- function(v, s) mns(paste0(num(v), " ± ", num(s)))  # v +/- s
sup  <- c("0"="⁰","1"="¹","2"="²","3"="³","4"="⁴",
          "5"="⁵","6"="⁶","7"="⁷","8"="⁸","9"="⁹")
sci  <- function(p){
  if(is.na(p)) return("—")
  if(p >= 5e-3) return(sprintf("%.3f", p))
  e <- floor(log10(p)); m <- p/10^e
  es <- paste0("⁻", paste(sup[strsplit(as.character(abs(e)),"")[[1]]], collapse=""))
  sprintf("%.1f×10%s", m, es)
}
## GxE test P: IID (shared) Wald test and FREE (heterogeneous) MV-Wald test
iidp <- function(s) if(is.list(s) && is.list(s$iid) && !is.null(s$iid$h2))
  Waldtest(s$iid$h2[2], s$iid$h2Covmat[2,2]) else NA
freep<- function(s){ if(!is.list(s) || !is.list(s$free) || is.null(s$free$h2)) return(NA)
  ng <- length(s$free$h2); het <- 2:(1+ng)
  tryCatch(MVWaldtest(s$free$sig2s[het], s$free$sig2Var[het,het]), error=function(e) NA) }

KCOLS <- c("liab_0.05"="Liability (K = 0.05)","liab_0.10"="Liability (K = 0.10)","liab_0.26"="Liability (K = 0.26)")

## Table 5
r  <- readRDS("reanalysis2022/cluster_gxemm/gxemm_cluster_results.rds")   # region, all 1,166 complete cases
fa <- readRDS("reanalysis2022/cluster_gxemm/gxemm_farmadj_results.rds")   # Model A = same region model, 1,120 subset
fm <- readRDS("reanalysis2022/cluster_gxemm/gxemm_farm_results.rds")      # farm, 1,120 (14 farms)
nfmt <- function(n) formatC(n, format = "d", big.mark = ",")
nR <- nfmt(1166); nS <- nfmt(fa$meta$N); stopifnot(fa$meta$N == fm$N)

regrow <- function(i, lab, obj = r, env = "Region", n = nR)       # per-region FREE h2 row
  c(env, n, lab, pm(obj$obs$free$h2[i], se(obj$obs$free$h2Covmat)[i]),
    sapply(names(KCOLS), function(k){ s<-obj[[k]]$free; pm(s$h2[i], se(s$h2Covmat)[i]) }))

prow <- function(env, n, obj, lab, fun, obscell=NULL)       # GxE p-value row
  c(env, n, lab, if(!is.null(obscell)) obscell else sci(fun(obj$obs)),
    sapply(names(KCOLS), function(k) if(is.null(obj[[k]])) "re-running" else sci(fun(obj[[k]]))))

## farm cells say "re-running" if that K is missing from the farm run
fsummary <- sapply(names(KCOLS), function(k)
  if(is.null(fm[[k]])) "re-running" else "non-identifiable (Supplementary Table 5)")

rows <- rbind(
  regrow(1,"Northeast — h² ± SE"),
  regrow(2,"Midwest — h² ± SE"),
  regrow(3,"South — h² ± SE"),
  prow("Region", nR, r, "GxE test — IID (shared) P", iidp),
  prow("Region", nR, r, "GxE test — FREE (heterog.) P", freep),
  prow("Region", nS, fa$A_region_adj, "GxE test — IID (shared) P", iidp),     # same model, 46 fewer deer
  prow("Region", nS, fa$A_region_adj, "GxE test — FREE (heterog.) P", freep),
  regrow(1,"Northeast — h² ± SE", fa$B_farm_adj, "Region (+ farm)", nS),     # region interaction, farm-adjusted
  regrow(2,"Midwest — h² ± SE",   fa$B_farm_adj, "Region (+ farm)", nS),
  regrow(3,"South — h² ± SE",     fa$B_farm_adj, "Region (+ farm)", nS),
  prow("Region (+ farm)", nS, fa$B_farm_adj, "GxE test — IID (shared) P", iidp),
  prow("Region (+ farm)", nS, fa$B_farm_adj, "GxE test — FREE (heterog.) P", freep),
  c("Farm", nS, "Per-farm h² ± SE (14 farms)","FREE did not converge", fsummary),
  prow("Farm", nS, fm, "GxE test — IID (shared) P", iidp),
  prow("Farm", nS, fm, "GxE test — FREE (heterog.) P", freep, obscell="—")
)
colnames(rows) <- c("Environment","n","Component","Observed", unname(KCOLS))
write.csv(rows, file.path(OUT,"table5_gxemm.csv"), row.names=FALSE)   # quoted: n has commas
cat("wrote", file.path(OUT,"table5_gxemm.csv"), "\n")

## Supplementary Table 5: per-farm heritability
cov <- read.csv("reanalysis2022/data/cov_1170.csv", stringsAsFactors=FALSE)
cov$cwd <- suppressWarnings(as.numeric(cov$binary)); cov$sexn <- suppressWarnings(as.numeric(cov$sex))
covf<- cov[!is.na(cov$cwd) & !is.na(cov$sexn) & !is.na(cov$age),]
farms <- fm$farms
sdf <- data.frame(
  Farm   = gsub("East","Northeast", farms),
  Region = gsub("[0-9]+","", gsub("East","Northeast", farms)),
  n      = sapply(farms, function(f) sum(covf$farm==f)),
  Prevalence = sapply(farms, function(f) sprintf("%.0f%%", 100*mean(covf$cwd[covf$farm==f]))),
  check.names=FALSE, stringsAsFactors=FALSE)
for(k in names(KCOLS)) if(is.list(fm[[k]]) && is.list(fm[[k]]$free)){
  s <- fm[[k]]$free
  sdf[[sprintf("Liability h² ± SE (K = %s)", sub("liab_","",k))]] <- mapply(pm, s$h2, se(s$h2Covmat))
}
sdf <- sdf[order(sdf$Region, as.integer(gsub("[^0-9]", "", sdf$Farm))),]   # numeric farm order
write.csv(sdf, file.path(OUT,"tableS5_gxemm_perfarm.csv"), row.names=FALSE, quote=FALSE)
cat("wrote", file.path(OUT,"tableS5_gxemm_perfarm.csv"), " (", nrow(sdf), "farms, K:",
    paste(sub("liab_","", names(KCOLS)[sapply(names(KCOLS), function(k) is.list(fm[[k]]))]), collapse=","), ")\n")
