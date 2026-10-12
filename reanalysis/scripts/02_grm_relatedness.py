#!/usr/bin/env python3
## Step 1b: relatedness from the GRM and family structure, plus PLINK IBD (PI_HAT)
## GRM off-diagonal is about 2 x kinship: 0.5 = 1st degree (parent-offspring, full sibs), 0.25 = 2nd degree, 0.125 = 3rd degree (first cousins)
## writes reanalysis/figures/fig4_grm_relatedness.png
import numpy as np, pandas as pd, gzip
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

## run from the main project folder
ROOT = "reanalysis"
S=f"{ROOT}/01_structure"; F=f"{ROOT}/figures"

## read in sample info and GRM
meta=pd.read_csv(f"{ROOT}/00_data/samples.tsv",sep="\t",dtype={"clean_id":str})
c2reg=dict(zip(meta.clean_id,meta.region)); c2bin=dict(zip(meta.clean_id,meta.binary))

ids=[l.split()[0] for l in open(f"{S}/grm.rel.id")]
n=len(ids)
G=np.loadtxt(gzip.open(f"{S}/grm.rel.gz","rt"))
assert G.shape==(n,n), G.shape
iu=np.triu_indices(n,k=1)
off=G[iu]
diag=np.diag(G)

print(f"GRM {n}x{n}; {len(off):,} unique off-diagonal pairs")
print("\n=== Diagonal (1+F inbreeding proxy) ===")
print(f"  mean={diag.mean():.3f}  sd={diag.std():.3f}  min={diag.min():.3f}  max={diag.max():.3f}")
print("\n=== Off-diagonal relatedness distribution ===")
for q in [0,1,5,25,50,75,95,99,99.9,100]:
    print(f"  {q:5.1f}th pct: {np.percentile(off,q): .4f}")
print(f"  mean={off.mean():.4f}  sd={off.std():.4f}")

## pair counts by GRM band
bands=[("duplicate/twin (>0.9)",0.9,9),("1st-degree (0.4-0.9)",0.4,0.9),
       ("2nd-degree (0.2-0.4)",0.2,0.4),("3rd-degree (0.1-0.2)",0.1,0.2),
       ("distant (0.05-0.1)",0.05,0.1)]
print("\n=== Related pair counts (GRM off-diagonal bands) ===")
for lab,lo,hi in bands:
    m=(off>=lo)&(off<hi); print(f"  {lab:26s}: {m.sum():6d} pairs")

## relatives per individual at GRM >= 0.2 (2nd degree or closer)
rel_thr=0.2
A=(G>=rel_thr).astype(int); np.fill_diagonal(A,0)
n_rel=A.sum(1)
print(f"\n=== Family structure: relatives per individual at GRM>={rel_thr} (>=2nd degree) ===")
print(f"  individuals with >=1 close relative: {(n_rel>=1).sum()} / {n} ({(n_rel>=1).mean()*100:.1f}%)")
print(f"  mean relatives/indiv={n_rel.mean():.2f}  max={n_rel.max()}")

## families = connected components of the GRM >= 0.2 graph (union-find)
parent=list(range(n))
def find(x):
    while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
    return x
def union(a,b):
    ra,rb=find(a),find(b)
    if ra!=rb: parent[ra]=rb
for a,b in zip(*np.where(np.triu(A,1))):
    union(a,b)
comp={}
for i in range(n): comp.setdefault(find(i),[]).append(i)
fams=[v for v in comp.values() if len(v)>=2]
fams.sort(key=len,reverse=True)
print(f"\n=== Families (connected components, GRM>={rel_thr}) ===")
print(f"  {len(fams)} families with >=2 members; largest sizes: {[len(f) for f in fams[:12]]}")
singl=sum(1 for v in comp.values() if len(v)==1)
print(f"  singletons (no close relative): {singl}")
## regions and CWD+ rate of the largest families
print("\n  Largest families: size | regions(count) | CWD+ rate")
for f in fams[:10]:
    regs=pd.Series([c2reg[ids[i]] for i in f]).value_counts().to_dict()
    br=np.mean([c2bin[ids[i]] for i in f])
    print(f"    size={len(f):3d}  regions={regs}  CWD+ rate={br:.2f}")

## mean relatedness within vs across regions
reg=np.array([c2reg[i] for i in ids])
same=reg[iu[0]]==reg[iu[1]]
print("\n=== Mean off-diagonal relatedness: within vs across region ===")
print(f"  within-region : mean={off[same].mean():.4f} (n={same.sum():,})")
print(f"  across-region : mean={off[~same].mean():.4f} (n={(~same).sum():,})")
for r in ["South","Midwest","East"]:
    mm=(reg[iu[0]]==r)&(reg[iu[1]]==r)
    print(f"    within {r:8s}: mean GRM={off[mm].mean():.4f}  (>0.2 pairs: {(off[mm]>0.2).sum()})")

## PLINK IBD (PI_HAT)
gen=pd.read_csv(f"{S}/ibd.genome",sep=r"\s+")
print("\n=== PLINK IBD PI_HAT relative-pair counts ===")
for lab,lo,hi in [("dup/MZ (>0.9)",0.9,9),("1st-deg (0.35-0.65)",0.35,0.65),
                  ("2nd-deg (0.15-0.35)",0.15,0.35),("3rd-deg (0.08-0.15)",0.08,0.15)]:
    m=(gen.PI_HAT>=lo)&(gen.PI_HAT<hi); print(f"  {lab:22s}: {m.sum():6d} pairs")
gen["reg1"]=gen.IID1.map(c2reg); gen["reg2"]=gen.IID2.map(c2reg)
hi1=gen[gen.PI_HAT>=0.35]
print(f"  1st-degree pairs within-region: {(hi1.reg1==hi1.reg2).sum()} / {len(hi1)}")

## plot GRM off-diagonal histogram
fig,ax=plt.subplots(1,2,figsize=(12,4.5))
ax[0].hist(off,bins=200,color="#444"); ax[0].set_yscale("log")
ax[0].axvline(0.2,ls="--",c="orange"); ax[0].axvline(0.4,ls="--",c="red")
ax[0].set_xlabel("GRM off-diagonal (~2*kinship)");ax[0].set_ylabel("pair count (log)")
ax[0].set_title("Genome-wide relatedness distribution (all 325,221 pairs)")
ax[1].hist(off[off>0.05],bins=120,color="#444")
ax[1].axvline(0.2,ls="--",c="orange",label="2nd-degree")
ax[1].axvline(0.4,ls="--",c="red",label="1st-degree")
ax[1].set_xlabel("GRM off-diagonal");ax[1].set_ylabel("pair count")
ax[1].set_title("Related pairs only (GRM>0.05)");ax[1].legend(frameon=False)
fig.tight_layout(); fig.savefig(f"{F}/fig4_grm_relatedness.png",dpi=150); plt.close(fig)
print("\nfigure written: fig4_grm_relatedness.png")
