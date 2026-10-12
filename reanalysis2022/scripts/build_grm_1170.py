#!/usr/bin/env python3
## VanRaden GRM (autosomes, MAF > 0.01) for the 2022 cohort, plus the covariate table cov_1170.csv read by every 2022 script
## (the analyses use the PLINK GRM from plink_2022.sh, not the K_1170.npy made here)
## writes reanalysis2022/data/K_1170.npy and cov_1170.csv
import os, sys, numpy as np, pandas as pd

## run from the main project folder, with the genotypes on stdin:
## unzip -p doi_10_5061_dryad_wh70rxwnt__v20220622.zip 1170_WTD_GWAAs_IncludingGxE.tsv | python3 reanalysis2022/scripts/build_grm_1170.py
D="doi_10_5061_dryad_wh70rxwnt__v20220622"
OUT="reanalysis2022/data"

## autosomal markers from the MarkerMap
mm=pd.read_csv(f"{D}/MarkerMap.tsv",sep="\t")
auto=set(mm.loc[mm.Chromosome.astype(str).isin([str(i) for i in range(1,30)]),"Variant"])

## read genotypes from stdin
hdr=sys.stdin.readline().rstrip("\n").split("\t")
meta_cols=hdr[:4]  # CMSID, GeneralRegion, Sex, Binary
auto_idx=[i for i,m in enumerate(hdr) if i>=4 and m in auto]
markers=[hdr[i] for i in auto_idx]
print(f"[grm] {len(auto_idx)} autosomal marker columns",flush=True)
ids=[]; regions=[]; sexes=[]; binaries=[]; rows=[]
for line in sys.stdin:
    f=line.rstrip("\n").split("\t")
    ids.append(f[0]); regions.append(f[1]); sexes.append(f[2]); binaries.append(f[3])
    g=[f[i] for i in auto_idx]
    rows.append(g)
    if len(ids)%200==0: print(f"[grm] {len(ids)} deer",flush=True)
n=len(ids); print(f"[grm] {n} deer x {len(markers)} SNPs -> encoding",flush=True)

## to int8 matrix, missing '?' as -1
G=np.empty((n,len(markers)),dtype=np.int8)
for i,r in enumerate(rows):
    G[i]=[(-1 if x=="?" else int(x)) for x in r]
del rows

## VanRaden GRM: drop MAF <= 0.01, mean-impute missing, standardise per marker
Gf=G.astype(np.float32); Gf[G<0]=np.nan
p=np.nanmean(Gf,0)/2.0
keep=(p>0.01)&(p<0.99)&np.isfinite(p)   # drop monomorphic/all-missing
Gf=Gf[:,keep]; p=p[keep]
sd=np.sqrt(2*p*(1-p))
# impute and standardise
inds=np.where(np.isnan(Gf)); Gf[inds]=np.take(2*p,inds[1])
Zs=((Gf-2*p)/sd).astype(np.float32)
M=Zs.shape[1]
K=(Zs@Zs.T)/M
print(f"[grm] GRM {K.shape} from {M} SNPs; diag mean={np.diag(K).mean():.3f}",flush=True)
np.save(f"{OUT}/K_1170.npy",K)

## covariate table, with age and farm from the metadata
md=pd.read_csv(f"{D}/metadata_1170_clean.csv").set_index("CMSID")
cov=pd.DataFrame({"CMSID":ids,"region":regions,"sex":sexes,"binary":binaries})
cov["age"]=cov.CMSID.map(md["Age_OG"]); cov["farm"]=cov.CMSID.map(md["Region_Farm"])
cov.to_csv(f"{OUT}/cov_1170.csv",index=False)
print(f"[grm] saved K_1170.npy and cov_1170.csv (n={n})",flush=True)
