#!/usr/bin/env python3
## VanRaden GRM (autosomes, MAF > 0.01) for the 807 deer of the 2020 cohort, same steps as build_grm_1170.py,
## plus the covariate table cov_807.csv (the analyses use the PLINK GRM from 00b_plink_2020.sh, not this GRM)
## writes K_807.npy, K_807.bin and cov_807.csv to reanalysis/01_structure
import os, numpy as np, pandas as pd

## run from the main project folder
D   = "doi_10_5061_dryad_xd2547dcw__v20200304"
OUT = "reanalysis/01_structure"
GENO = f"{D}/807WTD_DRYAD.csv"

## autosomal markers (chr 1-29) from the MarkerMap
mm = pd.read_csv(f"{D}/807WTD_MarkerMap.csv")
auto = set(mm.loc[mm.Chromosome.astype(str).isin([str(i) for i in range(1, 30)]), "Variant"])

## read genotypes
with open(GENO) as fh:
    hdr = fh.readline().rstrip("\n").split(",")
    meta = hdr[:6]                       # SampleID, GeneralRegion, Sex, Age_OG, CMS_CWDScore, CWD_Binary
    auto_idx = np.array([i for i, m in enumerate(hdr) if i >= 6 and m in auto])
    print(f"[grm807] {len(auto_idx)} autosomal marker columns of {len(hdr)-6} total", flush=True)
    ids = []; region = []; sex = []; age = []; score = []; binary = []; g3 = []
    for line in fh:
        f = line.rstrip("\n").split(",")
        ids.append(f[0]); region.append(f[1]); sex.append(f[2]); age.append(f[3]); score.append(f[4]); binary.append(f[5])
        g3.append(np.array(f, dtype="<U3")[auto_idx])
        if len(ids) % 200 == 0: print(f"[grm807] {len(ids)} deer", flush=True)

## genotypes are stored as nucleotides ('A_G', '?_?' missing), so count one allele per marker (0/1/2)
## the GRM does not depend on which allele is counted, here the larger allele code
n = len(ids); G3 = np.vstack(g3); del g3
print(f"[grm807] {n} deer x {G3.shape[1]} markers -> nucleotide decode", flush=True)
V  = G3.view("<U1").reshape(n, G3.shape[1], 3)          # 'A_G' as [A, _, G]
c0 = V[:, :, 0]; c2 = V[:, :, 2]; del G3, V
def code(a):                                             # A/C/G/T as 1/2/3/4, missing as 0
    out = np.zeros(a.shape, dtype=np.int8)
    for ch, c in (("A", 1), ("C", 2), ("G", 3), ("T", 4)): out[a == ch] = c
    return out
c0c = code(c0); c2c = code(c2); del c0, c2
miss = (c0c == 0)
refc = np.maximum(c0c, c2c); refc[miss] = 0             # larger allele code per cell
refc = refc.max(axis=0)                                 # counted allele per marker
G = (c0c == refc).astype(np.int8) + (c2c == refc).astype(np.int8)
G[miss] = -1; del c0c, c2c, miss

## VanRaden GRM, as in build_grm_1170.py
Gf = G.astype(np.float32); Gf[G < 0] = np.nan
p  = np.nanmean(Gf, 0) / 2.0
keep = (p > 0.01) & (p < 0.99) & np.isfinite(p)
Gf = Gf[:, keep]; p = p[keep]; sd = np.sqrt(2 * p * (1 - p))
inds = np.where(np.isnan(Gf)); Gf[inds] = np.take(2 * p, inds[1])
Zs = ((Gf - 2 * p) / sd).astype(np.float32); M = Zs.shape[1]
K = (Zs @ Zs.T) / M
print(f"[grm807] GRM {K.shape} from {M} SNPs (MAF>0.01); diag mean={np.diag(K).mean():.3f}", flush=True)

## save GRM and covariate table (cluster_gxemm/build_K_from_plink.R later overwrites K_807.bin with the PLINK GRM)
np.save(f"{OUT}/K_807.npy", K)
K.astype(np.float64).tofile(f"{OUT}/K_807.bin")         # raw float64 for R readBin
cov = pd.DataFrame({"SampleID": ids, "region": region, "sex": sex, "age": age, "score": score, "binary": binary})
cov.to_csv(f"{OUT}/cov_807.csv", index=False)
print(f"[grm807] saved K_807.npy, K_807.bin, cov_807.csv (n={n}, M={M})", flush=True)
