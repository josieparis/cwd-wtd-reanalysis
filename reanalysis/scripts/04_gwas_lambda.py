#!/usr/bin/env python3
## Step 2a/2b: GWAS of both traits (binary 0/1 and score 0/1/2, both fitted as linear like the authors' EMMAX) and lambda
## models: naive OLS, OLS + 10 genomic PCs and EMMAX (mixed model with the GRM), all intercept only
## writes reanalysis/01_structure/gwas_{binary,score}.tsv.gz and lambda_summary.json
import numpy as np, pandas as pd, gzip, json
from scipy import stats, optimize

## run from the main project folder
ROOT = "reanalysis"
S=f"{ROOT}/01_structure"
SCR = "reanalysis/work"   # plink --recode A files (00b_plink_2020.sh)
CHI2_MED=stats.chi2.ppf(0.5,1)   # 0.4549

## genomic control lambda
def lam(p):
    p=np.asarray(p); p=p[np.isfinite(p)&(p>0)]
    return np.median(stats.chi2.isf(p,1))/CHI2_MED

## read in genotypes (additive .raw)
print("loading genotypes ...",flush=True)
raw=pd.read_csv(f"{SCR}/deer_auto.raw",sep=r"\s+")
ids=raw["IID"].astype(str).values
snpcols=list(raw.columns[6:])
Z=raw[snpcols].to_numpy(np.float32); del raw
markers=np.array([c.rsplit("_",1)[0] for c in snpcols])
effA   =np.array([c.rsplit("_",1)[1] for c in snpcols])
# mean impute missing genotypes
cm=np.nanmean(Z,axis=0); nanpos=np.where(np.isnan(Z)); Z[nanpos]=np.take(cm,nanpos[1])
n,m=Z.shape; print(f"  {n} samples x {m} SNPs",flush=True)

## phenotypes, PCs and GRM in genotype order
meta=pd.read_csv(f"{ROOT}/00_data/samples.tsv",sep="\t",dtype={"clean_id":str}).set_index("clean_id")
meta=meta.loc[ids]
ybin=meta["binary"].to_numpy(float); yscr=meta["score"].to_numpy(float)
ev=pd.read_csv(f"{S}/pca.eigenvec",sep="\t"); ev.columns=["FID","IID"]+[f"PC{i}" for i in range(1,21)]
ev["IID"]=ev["IID"].astype(str); ev=ev.set_index("IID").loc[ids]
PC=ev[[f"PC{i}" for i in range(1,11)]].to_numpy(float)

gids=[l.split()[0] for l in open(f"{S}/grm.rel.id")]
G=np.loadtxt(gzip.open(f"{S}/grm.rel.gz","rt"))
gi={g:k for k,g in enumerate(gids)}; order=[gi[i] for i in ids]
K=G[np.ix_(order,order)]                      # GRM in genotype order
evals,U=np.linalg.eigh(K); evals=np.clip(evals,1e-6,None)

## association scans
def scan_ols(y,C):
    # OLS scan of each SNP with fixed covariates C (n x q, incl intercept)
    Q,_=np.linalg.qr(C)
    yr=y-Q@(Q.T@y)
    Zr=Z-Q@(Q.T@Z)
    Szz=np.einsum("ij,ij->j",Zr,Zr)
    Szy=Zr.T@yr
    b=Szy/Szz
    df=n-C.shape[1]-1
    rss=(yr@yr)-b*Szy
    se=np.sqrt((rss/df)/Szz)
    t=b/se; p=2*stats.t.sf(np.abs(t),df)
    return b,se,p

def reml_delta(y,X):
    # REML estimate of delta = Ve/Vg and Vg (fixed effects X + GRM), using the eigen decomposition of K
    Uty=U.T@y; UtX=U.T@X; nq=n-X.shape[1]
    def negll(ld):
        d=np.exp(ld); w=1.0/(evals+d)
        XtWX=(UtX*w[:,None]).T@UtX
        XtWy=(UtX*w[:,None]).T@Uty
        beta=np.linalg.solve(XtWX,XtWy)
        r=Uty-UtX@beta
        vg=(w*r*r).sum()/nq
        ll=-0.5*(nq*np.log(vg)+np.log(evals+d).sum()+np.linalg.slogdet(XtWX)[1]+nq)
        return -ll
    res=optimize.minimize_scalar(negll,bounds=(-10,10),method="bounded")
    d=np.exp(res.x); w=1.0/(evals+d)
    XtWX=(UtX*w[:,None]).T@UtX; beta=np.linalg.solve(XtWX,(UtX*w[:,None]).T@Uty)
    r=Uty-UtX@beta; vg=(w*r*r).sum()/nq
    return d,vg,w

def scan_emmax(y,X):
    # EMMAX: variance components fixed from the null model (X only), then a Wald test per SNP
    d,vg,w=reml_delta(y,X)
    Uty=U.T@y; UtX=U.T@X; UtZ=U.T@Z
    # weighted residuals on X
    XtWX=(UtX*w[:,None]).T@UtX; XtWX_inv=np.linalg.inv(XtWX)
    def wresid(M):    # M: n x k, rotated
        return M-UtX@(XtWX_inv@((UtX*w[:,None]).T@M))
    yr=wresid(Uty[:,None])[:,0]; Zr=wresid(UtZ)
    denom=np.einsum("ij,ij->j",Zr*w[:,None],Zr)
    num  =(Zr*w[:,None]).T@yr
    b=num/denom
    chi2=b*b*denom/vg
    p=stats.chi2.sf(chi2,1)
    se=np.sqrt(vg/denom)
    ph=1.0/(1.0+d)   # pseudo-heritability from the null
    return b,se,p,ph,d

## run scans for both traits and save
ones=np.ones((n,1))
Xpc=np.column_stack([ones,PC])
results={}
for trait,y in [("binary",ybin),("score",yscr)]:
    print(f"\n=== {trait} ===",flush=True)
    bo,seo,po=scan_ols(y,ones)
    bp,sep,pp=scan_ols(y,Xpc)
    be,see,pe,ph,delta=scan_emmax(y,ones)
    print(f"  lambda  naive={lam(po):.3f}  PC10={lam(pp):.3f}  EMMAX={lam(pe):.3f}  (EMMAX pseudo-h2={ph:.3f})")
    df=pd.DataFrame({"marker":markers,"effA":effA,
        "b_naive":bo,"p_naive":po,"b_pc":bp,"p_pc":pp,"b_emmax":be,"se_emmax":see,"p_emmax":pe})
    df.to_parquet(f"{S}/gwas_{trait}.parquet") if False else df.to_csv(f"{S}/gwas_{trait}.tsv.gz",sep="\t",index=False)
    results[trait]={"lambda":{"naive":lam(po),"pc10":lam(pp),"emmax":lam(pe)},
                    "emmax_pseudo_h2":float(ph),"emmax_delta":float(delta)}

json.dump(results,open(f"{S}/lambda_summary.json","w"),indent=2)
print("\nsaved gwas_{binary,score}.tsv.gz and lambda_summary.json")
print(json.dumps(results,indent=2))
