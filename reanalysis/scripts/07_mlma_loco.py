#!/usr/bin/env python3
## Step 2 sensitivity: MLMA-LOCO, the SNPs on each autosome are tested with EMMAX using a GRM from all the other
## autosomes, so a SNP is never in its own GRM; compares lambda and the authors' hits (PRNP and the rest) with single-GRM EMMAX
## writes reanalysis/01_structure/mlma_loco_{binary,score}.tsv.gz
import numpy as np, pandas as pd
from scipy import stats, optimize

## run from the main project folder
ROOT = "reanalysis"
S=f"{ROOT}/01_structure"
SCR = "reanalysis/work"   # plink --recode A files (00b_plink_2020.sh)
CHI2_MED=stats.chi2.ppf(0.5,1)
def lam(p): p=np.asarray(p); p=p[np.isfinite(p)&(p>0)]; return np.median(stats.chi2.isf(p,1))/CHI2_MED

## read in genotypes and mean impute missing
print("loading genotypes ...",flush=True)
raw=pd.read_csv(f"{SCR}/deer_auto.raw",sep=r"\s+")
ids=raw["IID"].astype(str).values
snpcols=list(raw.columns[6:])
Z=raw[snpcols].to_numpy(np.float32); del raw
markers=np.array([c.rsplit("_",1)[0] for c in snpcols])
cm=np.nanmean(Z,axis=0); nanpos=np.where(np.isnan(Z)); Z[nanpos]=np.take(cm,nanpos[1])
n,m=Z.shape

## chromosome of each SNP (from the bim)
bim=pd.read_csv(f"{S}/deer_auto.bim",sep="\t",header=None,names=["chr","marker","cmm","pos","a1","a2"])
c2chr=dict(zip(bim.marker,bim.chr)); chrom=np.array([c2chr[x] for x in markers])

## standardise genotypes (VanRaden): z = (x - 2p) / sqrt(2p(1-p)), p = mean/2
p=Z.mean(0)/2.0; sd=np.sqrt(2*p*(1-p)); sd[sd==0]=1
Zs=((Z-2*p)/sd).astype(np.float32); del Z

## per-chromosome GRM parts (Zs_c Zs_c^T, not normalised) and SNP counts
chrs=list(range(1,30))
Sc={}; Mc={}
Ksum=np.zeros((n,n))
for c in chrs:
    idx=np.where(chrom==c)[0]
    Zc=Zs[:,idx]
    Smat=(Zc@Zc.T).astype(np.float64)
    Sc[c]=Smat; Mc[c]=len(idx); Ksum+=Smat
Mtot=sum(Mc.values())
print(f"  {n} samples, {Mtot} autosomal SNPs across {len(chrs)} chromosomes",flush=True)

## phenotypes
meta=pd.read_csv(f"{ROOT}/00_data/samples.tsv",sep="\t",dtype={"clean_id":str}).set_index("clean_id").loc[ids]
ybin=meta["binary"].to_numpy(float); yscr=meta["score"].to_numpy(float)
ones=np.ones((n,1))

def reml_and_scan(y, K, snp_idx):
    # EMMAX under GRM K: null REML (intercept only), then a Wald test of the SNPs in snp_idx
    evals,U=np.linalg.eigh(K); evals=np.clip(evals,1e-6,None)
    Uty=U.T@y; UtX=U.T@ones; nq=n-1
    def negll(ld):
        d=np.exp(ld); w=1.0/(evals+d)
        XtWX=(UtX*w[:,None]).T@UtX; beta=np.linalg.solve(XtWX,(UtX*w[:,None]).T@Uty)
        r=Uty-UtX@beta; vg=(w*r*r).sum()/nq
        return 0.5*(nq*np.log(vg)+np.log(evals+d).sum()+np.linalg.slogdet(XtWX)[1]+nq)
    d=np.exp(optimize.minimize_scalar(negll,bounds=(-10,10),method="bounded").x)
    w=1.0/(evals+d)
    XtWX=(UtX*w[:,None]).T@UtX; beta=np.linalg.solve(XtWX,(UtX*w[:,None]).T@Uty)
    r=Uty-UtX@beta; vg=(w*r*r).sum()/nq
    XtWX_inv=np.linalg.inv(XtWX)
    UtZ=U.T@Zs[:,snp_idx]
    yr=(Uty-UtX@(XtWX_inv@((UtX*w[:,None]).T@Uty)))
    Zr=UtZ-UtX@(XtWX_inv@((UtX*w[:,None]).T@UtZ))
    denom=np.einsum("ij,ij->j",Zr*w[:,None],Zr); num=(Zr*w[:,None]).T@yr
    b=num/denom; chi2=b*b*denom/vg
    return b, stats.chi2.sf(chi2,1), 1.0/(1.0+d)

## LOCO scans: GRM without chromosome c, test the SNPs on c
out={}
for trait,y in [("binary",ybin),("score",yscr)]:
    pv=np.full(m,np.nan); bv=np.full(m,np.nan); h2s=[]
    for c in chrs:
        Kloco=(Ksum-Sc[c])/(Mtot-Mc[c])
        idx=np.where(chrom==c)[0]
        b,p_,ph=reml_and_scan(y,Kloco,idx)
        pv[idx]=p_; bv[idx]=b; h2s.append(ph)
    out[trait]=pd.DataFrame({"marker":markers,"chr":chrom,"b_loco":bv,"p_loco":pv})
    out[trait].to_csv(f"{S}/mlma_loco_{trait}.tsv.gz",sep="\t",index=False)
    print(f"\n{trait}: LOCO λ_GC = {lam(pv):.3f}   (mean LOCO pseudo-h2={np.mean(h2s):.3f})",flush=True)

## compare with single-GRM EMMAX for PRNP and the authors' hits
hits=pd.read_csv(f"{ROOT}/00_data/tableS1_hits.tsv",sep="\t")
for trait,models in [("binary",["Binary_NoCovariates"]),("score",["IntervalTissueScores_NoCovariates"])]:
    em=pd.read_csv(f"{S}/gwas_{trait}.tsv.gz",sep="\t")[["marker","p_naive","p_emmax"]]
    lo=out[trait][["marker","p_loco"]]
    h=hits[hits.trait_model.isin(models)].merge(em,on="marker",how="left").merge(lo,on="marker",how="left")
    for cc in ["p_naive","p_emmax","p_loco"]: h[cc]=pd.to_numeric(h[cc],errors="coerce")
    h=h.sort_values("p_emmax")
    fe=lambda x: f"{x:9.1e}" if pd.notna(x) else f"{'  NA':>9}"
    print(f"\n=== {trait}: reported hits — naive vs single-GRM EMMAX vs LOCO ===")
    print(f"{'gene':22s}{'marker':20s}{'p_naive':>10}{'p_EMMAX':>10}{'p_LOCO':>10}")
    for _,r in h.iterrows():
        print(f"{str(r['gene'])[:21]:22s}{str(r['marker'])[:19]:20s}{fe(r.p_naive)}{fe(r.p_emmax)}{fe(r.p_loco)}")
print("\nsaved mlma_loco_{binary,score}.tsv.gz")
