#!/usr/bin/env python3
## PLINK .tped/.tfam for the 2022 cohort (autosomes 1-29) from the Dryad genotype file on stdin, for the GRM
## genotypes are 0/1/2 ('?' missing), written with dummy alleles A/C; plink_2022.sh makes the .bed from these
## unzip -p doi_10_5061_dryad_wh70rxwnt__v20220622.zip 1170_WTD_GWAAs_IncludingGxE.tsv | python3 reanalysis2022/scripts/build_bed_2022.py OUTPREFIX
import os, sys, numpy as np, pandas as pd
OUT = sys.argv[1]
## run from the main project folder
MM  = "doi_10_5061_dryad_wh70rxwnt__v20220622/MarkerMap.tsv"

## autosomal markers
mm  = pd.read_csv(MM, sep="\t")
mm  = mm[mm.Chromosome.astype(str).isin([str(i) for i in range(1, 30)])]
pos = dict(zip(mm.Variant, zip(mm.Chromosome.astype(int), mm.Position.astype(int))))

## read genotypes from stdin
hdr = sys.stdin.readline().rstrip("\n").split("\t")           # CMSID, GeneralRegion, Sex, Binary, markers
mk_idx = [i for i in range(4, len(hdr)) if hdr[i] in pos]     # autosomal marker columns
markers = [hdr[i] for i in mk_idx]
sys.stderr.write(f"[bed] {len(markers)} autosomal marker columns\n")

ids = []; rows = []
for line in sys.stdin:
    f = line.rstrip("\n").split("\t"); ids.append(f[0])
    rows.append([f[i] for i in mk_idx])
n = len(ids); sys.stderr.write(f"[bed] {n} deer -> encode\n")
G = np.empty((n, len(markers)), dtype=np.int8)
for i, r in enumerate(rows):
    G[i] = [(-1 if x == "?" else int(x)) for x in r]
del rows

## write .tfam and .tped
gmap = np.array(["0 0", "A A", "A C", "C C"])                 # index G+1: missing, AA, AC, CC
with open(f"{OUT}.tfam", "w") as tf:
    for s in ids: tf.write(f"{s} {s} 0 0 0 -9\n")
with open(f"{OUT}.tped", "w") as tp:
    for j, snp in enumerate(markers):
        c, p = pos[snp]
        tp.write(f"{c}\t{snp}\t0\t{p}\t" + " ".join(gmap[G[:, j] + 1]) + "\n")
        if j % 20000 == 0: sys.stderr.write(f"[bed] wrote {j} markers\n")
sys.stderr.write(f"[bed] done: {OUT}.tped/.tfam ({n} deer x {len(markers)} markers)\n")
