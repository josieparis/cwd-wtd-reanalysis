#!/usr/bin/env python3
## Step 2b: inflation and attenuation; QQ plots (naive, PC, EMMAX), hit counts per model
## and the p-values of the authors' Table S1 hits under each model (PRNP vs the rest)
## writes reanalysis/figures/fig5_qq_naive_pc_emmax.png and fig6_manhattan_emmax_binary.png
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats

## run from the main project folder
ROOT = "reanalysis"
S=f"{ROOT}/01_structure"; F=f"{ROOT}/figures"

## genomic control lambda
CHI2_MED=stats.chi2.ppf(0.5,1)
def lam(p): p=p[np.isfinite(p)&(p>0)]; return np.median(stats.chi2.isf(p,1))/CHI2_MED

## thresholds: genome-wide, Bonferroni (120,113 SNPs) and suggestive
gw=5e-8; bonf=0.05/120113; sugg=1e-5
## read in GWAS results and SNP positions
gb=pd.read_csv(f"{S}/gwas_binary.tsv.gz",sep="\t")
gs=pd.read_csv(f"{S}/gwas_score.tsv.gz",sep="\t")
mmap=pd.read_csv(f"{ROOT}/00_data/deer.bim",sep="\t",header=None,names=["chr","marker","cm","pos","a1","a2"])
gb=gb.merge(mmap[["marker","chr","pos"]],on="marker",how="left")
gs=gs.merge(mmap[["marker","chr","pos"]],on="marker",how="left")

## SNP counts below the thresholds
print("=== SNP counts below thresholds (of 120,113 autosomal SNPs) ===")
print(f"{'trait/model':22s}{'p<1e-5':>9}{'p<Bonf(4.2e-7)':>16}{'p<5e-8':>9}   lambda")
for tag,df,col in [("binary NAIVE",gb,"p_naive"),("binary PC10",gb,"p_pc"),("binary EMMAX",gb,"p_emmax"),
                   ("score  NAIVE",gs,"p_naive"),("score  PC10",gs,"p_pc"),("score  EMMAX",gs,"p_emmax")]:
    p=df[col].values
    print(f"{tag:22s}{(p<sugg).sum():9d}{(p<bonf).sum():16d}{(p<gw).sum():9d}   {lam(p):.2f}")

## QQ plots
def qqdata(p):
    p=np.sort(p[np.isfinite(p)&(p>0)]); n=len(p)
    exp=-np.log10((np.arange(1,n+1)-0.5)/n); obs=-np.log10(p); return exp,obs
fig,ax=plt.subplots(1,2,figsize=(12,5.4))
for a,(df,ttl) in zip(ax,[(gb,"Binary CWD"),(gs,"CWD tissue score")]):
    for col,lab,c in [("p_naive","naive (no correction)","#d62728"),
                      ("p_pc","+10 genomic PCs","#ff7f0e"),
                      ("p_emmax","EMMAX (GRM mixed model)","#1f77b4")]:
        e,o=qqdata(df[col].values); a.plot(e,o,".",ms=2,c=c,label=f"{lab}  λ={lam(df[col].values):.2f}")
    m=max(e.max(),8); a.plot([0,m],[0,m],"k--",lw=.8)
    a.axhline(-np.log10(gw),ls=":",c="grey",lw=.8)
    a.set_xlabel("expected -log10(p)"); a.set_ylabel("observed -log10(p)")
    a.set_title(f"QQ: {ttl}"); a.legend(fontsize=8,frameon=False,markerscale=4)
fig.tight_layout(); fig.savefig(f"{F}/fig5_qq_naive_pc_emmax.png",dpi=150); plt.close(fig)

## Manhattan plot (EMMAX binary) with PRNP highlighted
d=gb.dropna(subset=["chr","pos"]).sort_values(["chr","pos"]).copy()
d=d[d["chr"]<=29]; d["i"]=range(len(d))
fig,a=plt.subplots(figsize=(13,4))
for ch,g in d.groupby("chr"):
    a.plot(g["i"],-np.log10(g["p_emmax"]),".",ms=1.5,c=("#666" if int(ch)%2 else "#bbb"))
prnp=d[d.marker=="G96S_G_to_A_at286"]
if len(prnp): a.plot(prnp["i"],-np.log10(prnp["p_emmax"]),"*",ms=15,c="red",label="PRNP G96S")
a.axhline(-np.log10(gw),ls="--",c="red",lw=.7); a.axhline(-np.log10(bonf),ls=":",c="blue",lw=.7)
a.set_ylabel("-log10(p) EMMAX"); a.set_xlabel("autosomal position (chr order)")
a.set_title("EMMAX binary GWAS (structure-corrected) — only PRNP survives"); a.legend(frameon=False)
fig.tight_layout(); fig.savefig(f"{F}/fig6_manhattan_emmax_binary.png",dpi=150); plt.close(fig)

## the authors' Table S1 hits under each model
hits=pd.read_csv(f"{ROOT}/00_data/tableS1_hits.tsv",sep="\t")
def track(df,models,label):
    h=hits[hits.trait_model.isin(models)].merge(
        df[["marker","p_naive","p_pc","p_emmax"]],on="marker",how="left")
    for c in ["p_reported","p_naive","p_pc","p_emmax"]:
        h[c]=pd.to_numeric(h[c],errors="coerce")
    h=h.sort_values("p_naive")
    fe=lambda x: f"{x:10.1e}" if pd.notna(x) else f"{'  NA(X/MT)':>10}"
    print(f"\n=== {label}: authors' reported hits across MY models ===")
    print(f"{'gene':22s}{'marker':20s}{'p_report':>10}{'p_naive':>10}{'p_PC':>10}{'p_EMMAX':>10}")
    for _,r in h.iterrows():
        g=str(r['gene'])[:21]
        print(f"{g:22s}{str(r['marker'])[:19]:20s}{fe(r['p_reported'])}{fe(r['p_naive'])}{fe(r['p_pc'])}{fe(r['p_emmax'])}")
    return h
track(gb,["Binary_NoCovariates"],"BINARY (no-cov EMMAX reported)")
track(gs,["IntervalTissueScores_NoCovariates"],"SCORE (no-cov EMMAX reported)")
print("\nfigures: fig5_qq_naive_pc_emmax.png, fig6_manhattan_emmax_binary.png")
