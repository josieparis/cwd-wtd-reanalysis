#!/usr/bin/env python3
## Step 2c: winner's curse on PRNP codon 96 (G96S), full sample (n = 807) vs many smaller "discovery" subsamples
## effects and R2 that pass a significance threshold in small samples are biased upwards
## writes reanalysis/01_structure/winners_curse.tsv, figures/fig7_winners_curse.png and fig8_winners_curse_scatter.png
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats

## run from the main project folder
ROOT = "reanalysis"
S=f"{ROOT}/01_structure"; F=f"{ROOT}/figures"
SCR = "reanalysis/work"   # plink --recode A files (00b_plink_2020.sh)
rng=np.random.default_rng(20200304)

## PRNP genotypes and phenotypes
raw=pd.read_csv(f"{SCR}/prnp.raw",sep=r"\s+").set_index("IID")
meta=pd.read_csv(f"{ROOT}/00_data/samples.tsv",sep="\t",dtype={"clean_id":str}).set_index("clean_id")
ids=[i for i in raw.index if i in meta.index]
g=raw.loc[ids,"G96S_G_to_A_at286_A"].to_numpy(float)     # dosage of the protective A (S96) allele
ok=~np.isnan(g); g=g[ok]
ybin=meta.loc[ids,"binary"].to_numpy(float)[ok]
yscr=meta.loc[ids,"score"].to_numpy(float)[ok]
N=len(g); print(f"PRNP G96S genotyped in {N} animals; MAF(A)={g.sum()/(2*N):.3f}")
print("genotype counts:",dict(zip(*np.unique(g,return_counts=True))))

def lin(y,x):
    # simple linear regression of y on x: slope, p, R2
    x=x-x.mean(); y=y-y.mean(); sxx=(x*x).sum()
    b=(x*y).sum()/sxx
    n=len(y); resid=y-b*x; s2=(resid*resid).sum()/(n-2)
    se=np.sqrt(s2/sxx); t=b/se; p=2*stats.t.sf(abs(t),n-2)
    r2=(np.corrcoef(x,y)[0,1])**2
    return b,p,r2

## full-sample effects
for tag,y in [("BINARY",ybin),("SCORE",yscr)]:
    b,p,r2=lin(y,g)
    lo=stats.logistic.cdf  # not used
    print(f"\nFull-sample {tag}: beta(per A allele)={b:.4f}  p={p:.2e}  R^2(var explained)={r2:.4f}")
## logistic OR for binary
import math
X=np.column_stack([np.ones(N),g]);
# Newton-Raphson logistic fit
beta=np.zeros(2)
for _ in range(50):
    eta=X@beta; mu=1/(1+np.exp(-eta)); Wv=mu*(1-mu)
    z=X.T@(ybin-mu); H=X.T@(X*Wv[:,None])
    beta+=np.linalg.solve(H,z)
se_or=np.sqrt(np.linalg.inv(X.T@(X*(1/(1+np.exp(-X@beta))*(1-1/(1+np.exp(-X@beta))))[:,None]))[1,1])
print(f"Full-sample BINARY logistic: OR per A allele = {np.exp(beta[1]):.3f} "
      f"(95%CI {np.exp(beta[1]-1.96*se_or):.3f}-{np.exp(beta[1]+1.96*se_or):.3f})")

## subsample smaller "discovery" studies (3000 draws per N)
Ngrid=[40,60,80,120,160,220,300,400,550,700]
B=3000
gw=5e-8; nom=0.05
def subsample(y,yname):
    full_b,full_p,full_r2=lin(y,g)
    rows=[]
    for n in Ngrid:
        bs=np.empty(B); ps=np.empty(B); r2s=np.empty(B)
        for k in range(B):
            idx=rng.choice(N,n,replace=False)
            bs[k],ps[k],r2s[k]=lin(y[idx],g[idx])
        for thr,lab in [(nom,"p<0.05"),(gw,"p<5e-8")]:
            sig=ps<thr
            rows.append(dict(y=yname,N=n,thr=lab,power=sig.mean(),
                mean_b_all=np.abs(bs).mean(), mean_b_sig=(np.abs(bs[sig]).mean() if sig.any() else np.nan),
                mean_r2_all=r2s.mean(), mean_r2_sig=(r2s[sig].mean() if sig.any() else np.nan),
                full_r2=full_r2, full_b=abs(full_b)))
    return pd.DataFrame(rows)

res=pd.concat([subsample(ybin,"binary"),subsample(yscr,"score")],ignore_index=True)
res.to_csv(f"{S}/winners_curse.tsv",sep="\t",index=False)

## R2 in the discovered subsamples vs the full sample
print("\n=== Winner's curse: R^2 estimated in smaller 'discovery' samples (conditional on significance) ===")
for yname in ["binary","score"]:
    fr=res[(res.y==yname)&(res.thr=="p<5e-8")]
    full=fr.full_r2.iloc[0]
    print(f"\n {yname}: full-sample (n={N}) R^2 = {full:.4f}")
    print(f"   {'N':>5}{'power':>8}{'E[R2|discovered]':>18}{'inflation x':>12}")
    for _,r in fr.iterrows():
        infl=r.mean_r2_sig/full if r.mean_r2_sig==r.mean_r2_sig else np.nan
        print(f"   {int(r.N):5d}{r.power:8.2f}{r.mean_r2_sig:18.4f}{infl:12.2f}")

## plot
fig,ax=plt.subplots(1,2,figsize=(13,5))
for yname,c in [("binary","#1f77b4"),("score","#d62728")]:
    fr=res[(res.y==yname)&(res.thr=="p<5e-8")]
    full=fr.full_r2.iloc[0]
    ax[0].plot(fr.N,fr.mean_r2_sig,"o-",c=c,label=f"{yname}: E[R²|discovered]")
    ax[0].axhline(full,ls="--",c=c,alpha=.6,label=f"{yname}: full-sample R²={full:.3f}")
ax[0].set_xlabel("discovery sample size N"); ax[0].set_ylabel("PRNP variance explained (R²)")
ax[0].set_title("Winner's curse: variance explained shrinks as N grows\n(estimates conditional on genome-wide 'discovery')")
ax[0].legend(fontsize=7,frameon=False)
## power curve
for yname,c in [("binary","#1f77b4"),("score","#d62728")]:
    fr=res[(res.y==yname)&(res.thr=="p<5e-8")]
    ax[1].plot(fr.N,fr.power,"o-",c=c,label=f"{yname} (p<5e-8)")
ax[1].set_xlabel("discovery sample size N"); ax[1].set_ylabel("P(PRNP reaches genome-wide sig.)")
ax[1].set_title("Power to 'discover' PRNP by sample size"); ax[1].legend(frameon=False)
fig.tight_layout(); fig.savefig(f"{F}/fig7_winners_curse.png",dpi=150); plt.close(fig)

## effect vs significance at a small N (n = 80)
n0=80; B2=1500; bb=np.empty(B2); pp=np.empty(B2)
for k in range(B2):
    idx=rng.choice(N,n0,replace=False); bb[k],pp[k],_=lin(ybin[idx],g[idx])
fig,a=plt.subplots(figsize=(6.5,5))
sig=pp<gw
a.scatter(np.abs(bb[~sig]),-np.log10(pp[~sig]),s=8,c="#bbb",label="not discovered")
a.scatter(np.abs(bb[sig]),-np.log10(pp[sig]),s=10,c="#d62728",label="'discovered' (p<5e-8)")
a.axhline(-np.log10(gw),ls="--",c="k",lw=.7)
a.axvline(abs(lin(ybin,g)[0]),ls="--",c="blue",lw=1,label=f"true |beta| (n={N})")
a.set_xlabel("|estimated PRNP effect| per A allele"); a.set_ylabel("-log10(p)")
a.set_title(f"At N={n0}: only large-effect draws are 'discovered'\n(selection on significance = inflated effects)")
a.legend(fontsize=8,frameon=False)
fig.tight_layout(); fig.savefig(f"{F}/fig8_winners_curse_scatter.png",dpi=150); plt.close(fig)
print("\nfigures: fig7_winners_curse.png, fig8_winners_curse_scatter.png")
