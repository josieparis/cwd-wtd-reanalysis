## genotype concordance (IBS allele sharing) between every 2020 deer and every 2022 deer on shared markers,
## and how many 2020 deer have a 2022 match >= 0.90 (data for Supplementary Fig 1)
## run from the main project folder; writes hist_pairwise.csv, max_per_animal.csv and summary.txt to SCR
import numpy as np, zipfile, io, time
t0=time.time()
RAW=f"{__import__('os').environ['SCR']}/geno2020.raw"   # SCR is the folder with the plink .raw file
ZIP="doi_10_5061_dryad_wh70rxwnt__v20220622.zip"
MEMBER="1170_WTD_GWAAs_IncludingGxE.tsv"

## 2020 genotypes (A1 dosages from the plink .raw)
import pandas as pd
raw=pd.read_csv(RAW, sep=r"\s+")
mk20=[c.rsplit('_',1)[0] for c in raw.columns[6:]]               # Affx ids
d20=raw.iloc[:,6:].to_numpy(dtype=float)                         # 807 x M, NaN = missing
d20=np.where(np.isnan(d20), -1, d20).astype(np.int8)
print(f"[2020] {d20.shape} markers={len(mk20)}  ({time.time()-t0:.0f}s)")

want=set(mk20); col20={m:i for i,m in enumerate(mk20)}

## 2022 genotypes, streamed from the Dryad zip (which allele is counted is unknown)
zf=zipfile.ZipFile(ZIP); fh=io.TextIOWrapper(zf.open(MEMBER),encoding="utf-8",newline="")
head=fh.readline().rstrip("\n").split("\t")
mk22=head[4:]
keepidx=[i for i,m in enumerate(mk22) if m in want]              # markers also in 2020
mk22k=[mk22[i] for i in keepidx]
absidx=[i+4 for i in keepidx]                                    # column index after the 4 metadata columns
rows=[]; region=[]; biny=[]
for line in fh:
    p=line.rstrip("\n").split("\t")
    region.append(p[1]); biny.append(p[3])
    rows.append([p[j] for j in absidx])
zf.close()
def enc(x): return -1 if x=='?' else int(x)
d22=np.array([[enc(v) for v in r] for r in rows], dtype=np.int8) # 2022 deer x K
print(f"[2022] {d22.shape} kept={len(mk22k)}  ({time.time()-t0:.0f}s)")

## align to the common markers
common=[m for m in mk22k if m in col20]
i20=[col20[m] for m in common]
c22={m:j for j,m in enumerate(mk22k)}; j22=[c22[m] for m in common]
A=d20[:,i20].astype(np.int16)      # 807 x K
B=d22[:,j22].astype(np.int16)      # 2022 deer x K
K=A.shape[1]; print(f"[common] {K} markers")

## allele frequencies, and flip the 2022 dosages where the other allele is counted
def freq(M):
    m=M>=0; s=np.where(m,M,0).sum(0); n=m.sum(0)*2
    return s/np.maximum(n,1), m.mean(0)
p20,call20=freq(A); p22,call22=freq(B)
flip = np.abs(p20-p22) > np.abs(p20-(1-p22))
Bf=B.copy(); mB=Bf>=0
Bf[:,flip]=np.where(mB[:,flip], 2-Bf[:,flip], -1)
## keep markers with MAF > 0.02, frequency outside 0.42-0.58 (clear orientation) and call rate > 0.95
keep = (p20>0.02)&(p20<0.98)&((p20<0.42)|(p20>0.58))&(call20>0.95)&(call22>0.95)
A=A[:,keep]; Bf=Bf[:,keep]; K=int(keep.sum())
print(f"[filter] kept {K} markers (dropped {int((~keep).sum())})")

## concordance from one-hot genotypes (allele matches over comparable markers)
def onehot(M):
    return [(M==g).astype(np.float32) for g in (0,1,2)], (M>=0).astype(np.float32)
(A0,A1,A2),MA=onehot(A); (B0,B1,B2),MB=onehot(Bf)
comp =MA@MB.T
D = (A0@B1.T + A1@B0.T + A1@B2.T + A2@B1.T) + 2.0*(A0@B2.T + A2@B0.T)   # sum of |a-b| over comparable markers
conc = 1.0 - D/(2.0*np.maximum(comp,1))    # IBS allele sharing, as plink DST
print(f"[conc] {conc.shape}  ({time.time()-t0:.0f}s)")

## summaries
flat=conc.ravel()
mx=conc.max(1)                              # best 2022 match for each 2020 deer
THR=0.90
nmatch=int((mx>=THR).sum())
print(f"2020 animals with a 2022 match >= {THR}: {nmatch}/807 ({100*nmatch/807:.1f}%)")
print(f"median conc of matched pairs (>= {THR}): {np.median(flat[flat>=THR]):.3f}")
print(f"median conc of other pairs   (<  {THR}): {np.median(flat[flat<THR]):.3f}")
sm=np.sort(mx)
print("per-animal max around the gap (sorted[20:34]):", np.round(sm[20:34],3))

## write figure data
edges=np.arange(0.40,1.0001,0.005)
h,_=np.histogram(flat,bins=edges); mids=(edges[:-1]+edges[1:])/2
np.savetxt(f"{__import__('os').environ['SCR']}/hist_pairwise.csv",
           np.c_[mids,h], delimiter=",", header="mid,count", comments="", fmt=["%.4f","%d"])
np.savetxt(f"{__import__('os').environ['SCR']}/max_per_animal.csv",
           mx, delimiter=",", header="maxconc", comments="", fmt="%.5f")
with open(f"{__import__('os').environ['SCR']}/summary.txt","w") as o:
    o.write(f"markers_used={K}\nnmatch_ge_0.90={nmatch}\npct={100*nmatch/807:.1f}\n"
            f"median_match={np.median(flat[flat>=THR]):.3f}\nmedian_other={np.median(flat[flat<THR]):.3f}\n"
            f"npairs={flat.size}\n")
print(f"done ({time.time()-t0:.0f}s)")
