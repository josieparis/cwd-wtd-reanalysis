#!/usr/bin/env python3
## GEMMA input files for the 2022 cohort (BIMBAM genotypes, phenotype, covariates) from the genotype file on stdin
## GEMMA makes its own centred GRM (-gk), a check on the rrBLUP/GCTA h2 that use the VanRaden GRM
## writes reanalysis2022/cluster_h2/gemma/{geno.bimbam, pheno.txt, covar_intercept.txt, covar_region.txt, covar_farm.txt}
import sys, csv, os
import numpy as np

## run from the main project folder:
## unzip -p doi_10_5061_dryad_wh70rxwnt__v20220622.zip 1170_WTD_GWAAs_IncludingGxE.tsv | python3 reanalysis2022/scripts/build_gemma_inputs.py
D   = "doi_10_5061_dryad_wh70rxwnt__v20220622"
OUT = "reanalysis2022/cluster_h2/gemma"
os.makedirs(OUT, exist_ok=True)

## autosomal markers and metadata (age, farm)
auto = set(r["Variant"] for r in csv.DictReader(open(f"{D}/MarkerMap.tsv"), delimiter="\t")
           if str(r["Chromosome"]) in {str(i) for i in range(1, 30)})
md = {r["CMSID"]: r for r in csv.DictReader(open(f"{D}/metadata_1170_clean.csv"))}

## read genotypes from stdin
hdr = sys.stdin.readline().rstrip("\n").split("\t")
idx = [i for i, m in enumerate(hdr) if i >= 4 and m in auto]
markers = [hdr[i] for i in idx]
ids = []; regs = []; sexes = []; bins = []; rows = []
for line in sys.stdin:
    f = line.rstrip("\n").split("\t")
    ids.append(f[0]); regs.append(f[1]); sexes.append(f[2]); bins.append(f[3])
    rows.append([f[i] for i in idx])
    if len(ids) % 200 == 0:
        print(f"[prep] {len(ids)} deer", file=sys.stderr, flush=True)

## genotype matrix, keep SNPs with MAF > 0.01
n = len(ids); M = len(markers)
G = np.empty((n, M), dtype=np.float32)
for i, r in enumerate(rows):
    G[i] = [np.nan if x == "?" else float(x) for x in r]
del rows
p = np.nanmean(G, 0) / 2.0
keep = np.where((p > 0.01) & (p < 0.99) & np.isfinite(p))[0]
G = G[:, keep]; markers = [markers[k] for k in keep]; M = len(markers)
print(f"[prep] {n} deer x {M} polymorphic autosomal SNPs", file=sys.stderr, flush=True)

## genotypes
with open(f"{OUT}/geno.bimbam", "w") as g:               # rs, alleleA, alleleB, dosages
    for j in range(M):
        vals = ",".join("NA" if np.isnan(v) else f"{v:.0f}" for v in G[:, j])
        g.write(f"{markers[j]},A,B,{vals}\n")

## phenotype (NA if missing)
with open(f"{OUT}/pheno.txt", "w") as fp:
    for b in bins:
        fp.write((b if b in ("0", "1") else "NA") + "\n")

## covariates: intercept only; + sex, age and region; + sex, age and farm
regset = ["East", "Midwest", "South"]                    # East = reference level
farms  = sorted({md[i]["Region_Farm"] for i in ids if i in md and md[i]["Region_Farm"]})
sx = lambda s: s if s in ("0", "1") else "NA"
with open(f"{OUT}/covar_intercept.txt", "w") as c0, \
     open(f"{OUT}/covar_region.txt", "w") as cr, \
     open(f"{OUT}/covar_farm.txt", "w") as cf:
    for k, i in enumerate(ids):
        age = md.get(i, {}).get("Age_OG", "NA"); fm = md.get(i, {}).get("Region_Farm", "")
        c0.write("1\n")
        rd = "\t".join("1" if regs[k] == r else "0" for r in regset[1:])
        cr.write(f"1\t{sx(sexes[k])}\t{age}\t{rd}\n")
        fd = "\t".join("1" if fm == ff else "0" for ff in farms[1:])
        cf.write(f"1\t{sx(sexes[k])}\t{age}\t{fd}\n")

print(f"[prep] wrote geno.bimbam ({M} SNPs), pheno.txt, 3 covar files "
      f"(region: 2 dummies; farm: {len(farms)-1} dummies for {len(farms)} farms)", file=sys.stderr)
