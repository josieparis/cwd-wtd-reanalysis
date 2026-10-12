#!/usr/bin/env python3
## per-SNP tables for the GWAS calibration figure (fig2_calibration.R), traits fitted linearly as in EMMAX
## 2022 naive and 10-PC scans, merged with the GCTA --mlma p-values (writes reanalysis2022/data/gwas2022_binary.tsv.gz)
## 2020 scans without and with the CWD score as a covariate (writes reanalysis/01_structure/collider_2020.tsv.gz)
import numpy as np, pandas as pd
from scipy import stats

## run from the main project folder; WORK has the PLINK --recode A files deer2022.raw and deer_auto.raw
BASE = "."
WORK = "work"
CHI2_MED = stats.chi2.ppf(0.5, 1)

## genomic control lambda
def lam(p):
    p = np.asarray(p); p = p[np.isfinite(p) & (p > 0)]
    return np.median(stats.chi2.isf(p, 1)) / CHI2_MED

## read a PLINK .raw file, mean-impute missing genotypes
def load_raw(path):
    df = pd.read_csv(path, sep=r"\s+")
    iid = df["IID"].astype(str).values
    markers = np.array([c.rsplit("_", 1)[0] for c in df.columns[6:]])
    Z = df.iloc[:, 6:].to_numpy(np.float32); del df
    cm = np.nanmean(Z, 0); nan = np.where(np.isnan(Z)); Z[nan] = np.take(cm, nan[1])
    return iid, markers, Z

## OLS of y on each SNP, adjusting for the covariates in C (incl. intercept)
def scan(y, Z, C):
    n = len(y); Q, _ = np.linalg.qr(C)
    yr = y - Q @ (Q.T @ y); Zr = Z - Q @ (Q.T @ Z)
    Szz = np.einsum("ij,ij->j", Zr, Zr); Szy = Zr.T @ yr
    b = Szy / Szz; dof = n - C.shape[1] - 1
    rss = (yr @ yr) - b * Szy; se = np.sqrt((rss / dof) / Szz)
    t = b / se; p = 2 * stats.t.sf(np.abs(t), dof)
    return b, p

## 2022: naive and 10-PC scans
print("[2022] loading genotypes ...", flush=True)
iid, mk, Z = load_raw(f"{WORK}/deer2022.raw")
cov = pd.read_csv(f"{BASE}/reanalysis2022/data/cov_1170.csv").set_index("CMSID")
pe = pd.read_csv(f"{BASE}/reanalysis2022/01_structure/pca2022.eigenvec", sep=r"\s+", header=None)
pe = pe.set_index(1)[list(range(2, 12))]; pe.columns = [f"PC{i}" for i in range(1, 11)]
y = pd.to_numeric(cov.loc[iid, "binary"], errors="coerce").to_numpy(float)
PC = pe.loc[iid].to_numpy(float)
ok = np.isfinite(y)                                        # drop the deer with no CWD status
y, Zc, PCc = y[ok], Z[ok], PC[ok]
one = np.ones((ok.sum(), 1))
_, p_naive = scan(y, Zc, one)
_, p_pc    = scan(y, Zc, np.hstack([one, PCc]))
## merge with the GCTA --mlma p-values (PLINK GRM) and save
mlma = pd.read_csv(f"{BASE}/reanalysis2022/mlma_results/mlma_bin_plink.mlma", sep="\t").set_index("SNP")
g = pd.DataFrame({"marker": mk, "p_naive": p_naive, "p_pc": p_pc})
g["chr"] = mlma.reindex(g.marker)["Chr"].values; g["bp"] = mlma.reindex(g.marker)["bp"].values
g["p_grm"] = mlma.reindex(g.marker)["p"].values
g.to_csv(f"{BASE}/reanalysis2022/data/gwas2022_binary.tsv.gz", sep="\t", index=False)
prnp = g[g.marker == "G96S_G_to_A_at286"].iloc[0]
print(f"[2022] lambda naive={lam(p_naive):.2f} pc={lam(p_pc):.2f} grm={lam(g.p_grm):.2f}", flush=True)
print(f"[2022] SNPs<5e-5: naive={int((p_naive<5e-5).sum())} pc={int((p_pc<5e-5).sum())} grm={int((g.p_grm<5e-5).sum())}", flush=True)
print(f"[2022] PRNP p: naive={prnp.p_naive:.2e} pc={prnp.p_pc:.2e} grm={prnp.p_grm:.2e}", flush=True)

## 2020 collider scan: binary ~ sex + age + region, without and with the CWD score
print("[2020] loading genotypes ...", flush=True)
iid20, mk20, Z20 = load_raw(f"{WORK}/deer_auto.raw")
sm = pd.read_csv(f"{BASE}/reanalysis/00_data/samples.tsv", sep="\t", dtype={"clean_id": str}).set_index("clean_id").loc[iid20]
bim = pd.read_csv(f"{BASE}/reanalysis/00_data/deer.bim", sep="\t", header=None,
                  names=["chr", "marker", "cm", "pos", "a1", "a2"]).set_index("marker")
y20 = sm["binary"].to_numpy(float)
sex = sm["sex_raw"].to_numpy(float); age = sm["age"].to_numpy(float); score = sm["score"].to_numpy(float)
reg = pd.get_dummies(sm["region"], drop_first=True).to_numpy(float)   # 2 region dummies
one20 = np.ones((len(y20), 1))
C_base  = np.hstack([one20, sex[:, None], age[:, None], reg])
C_score = np.hstack([C_base, score[:, None]])
_, p_base  = scan(y20, Z20, C_base)
_, p_score = scan(y20, Z20, C_score)
## save
c = pd.DataFrame({"marker": mk20, "p_base": p_base, "p_score": p_score})
c["chr"] = bim.reindex(c.marker)["chr"].values
c.to_csv(f"{BASE}/reanalysis/01_structure/collider_2020.tsv.gz", sep="\t", index=False)
pr = c[c.marker == "G96S_G_to_A_at286"].iloc[0]
top_score = c.sort_values("p_score").iloc[0]
print(f"[2020 collider] PRNP: p_base={pr.p_base:.2e} -> p_score={pr.p_score:.3f}", flush=True)
print(f"[2020 collider] top hit WITH score: {top_score.marker} chr{int(top_score.chr)} p_score={top_score.p_score:.2e} (its p_base={top_score.p_base:.3f})", flush=True)
print("saved gwas2022_binary.tsv.gz and collider_2020.tsv.gz", flush=True)
