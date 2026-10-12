#!/usr/bin/env python3
## Step 0: build the 2020 dataset from the Dryad files (807 deer x 124,020 SNPs)
## writes reanalysis/00_data/samples.tsv (phenotypes + covariates), deer.map, dup_pairs.tsv (AF31 pairs) and reanalysis/work/deer.ped
## the .ped phenotype is left as -9, the binary and score phenotypes are taken from samples.tsv downstream
import csv, sys, os

## run from the main project folder
DATA = "doi_10_5061_dryad_xd2547dcw__v20200304"
GENO = os.path.join(DATA, "807WTD_DRYAD.csv")
MMAP = os.path.join(DATA, "807WTD_MarkerMap.csv")
AF31 = os.path.join(DATA, "AF31_IBSDistanceDuplicatesControlsTwins.csv")

ROOT = "reanalysis"
SCRATCH = "reanalysis/work"   # large intermediates
os.makedirs(SCRATCH, exist_ok=True)
OUT_META = os.path.join(ROOT, "00_data", "samples.tsv")
OUT_DUP  = os.path.join(ROOT, "00_data", "dup_pairs.tsv")
OUT_PED  = os.path.join(SCRATCH, "deer.ped")   # large, keep off Dropbox
OUT_MAP  = os.path.join(ROOT, "00_data", "deer.map")

## recode chromosomes: 1-29 stay, X is 30, everything else (MT, scaffolds, unplaced) is 0
def recode_chrom(c):
    c = c.strip()
    if c.isdigit():
        n = int(c)
        return str(n) if 1 <= n <= 29 else "0"
    if c == "X":
        return "30"
    return "0"

## 1. marker map {variant: (chrom, pos)}, matched by Affx id
print("[1/4] reading marker map ...", flush=True)
mm = {}
with open(MMAP, newline="") as fh:
    r = csv.reader(fh)
    next(r)
    for row in r:
        if len(row) < 3:
            continue
        var, chrom, pos = row[0], row[1], row[2]
        mm[var] = (recode_chrom(chrom), pos if pos.strip() else "0")
print(f"      {len(mm):,} markers in map", flush=True)

## 2. genotype header (marker order), write .map
print("[2/4] reading genotype header + writing .map ...", flush=True)
fgeno = open(GENO, newline="")
reader = csv.reader(fgeno)
header = next(reader)
meta_cols = header[:6]
markers   = header[6:]
n_markers = len(markers)
assert meta_cols == ["SampleID","GeneralRegion","Sex","Age_OG","CMS_CWDScore","CWD_Binary"], meta_cols

n_in_map = 0
n_missing_map = 0
with open(OUT_MAP, "w") as fmap:
    for m in markers:
        if m in mm:
            chrom, pos = mm[m]
            n_in_map += 1
        else:
            chrom, pos = "0", "0"
            n_missing_map += 1
        fmap.write(f"{chrom}\t{m}\t0\t{pos}\n")
print(f"      {n_markers:,} SNP columns; {n_in_map:,} matched map, "
      f"{n_missing_map:,} unmatched (chrom 0)", flush=True)

## 3. stream genotypes (coded like A_G, missing ?_?) to .ped and samples.tsv
print("[3/4] streaming genotypes -> .ped + samples.tsv ...", flush=True)
GT = {"A_A":"A A","C_C":"C C","G_G":"G G","T_T":"T T",
      "A_C":"A C","A_G":"A G","A_T":"A T","C_G":"C G","C_T":"C T","G_T":"G T",
      "C_A":"C A","G_A":"G A","T_A":"T A","G_C":"G C","T_C":"T C","T_G":"T G",
      "?_?":"0 0"}

fmeta = open(OUT_META, "w")
fmeta.write("clean_id\tsample_id\tregion\tsex_raw\tage\tscore\tbinary\n")
fped = open(OUT_PED, "w")

n = 0
bad_tokens = {}
for row in reader:
    n += 1
    cid = f"S{n:04d}"
    sample_id, region, sex, age, score, binary = row[:6]
    fmeta.write(f"{cid}\t{sample_id}\t{region}\t{sex}\t{age}\t{score}\t{binary}\n")
    # .ped: FID IID PAT MAT SEX PHENO genotypes
    out = [cid, cid, "0", "0", "0", "-9"]
    geno = row[6:]
    for g in geno:
        v = GT.get(g)
        if v is None:
            bad_tokens[g] = bad_tokens.get(g, 0) + 1
            v = "0 0"
        out.append(v)
    fped.write(" ".join(out))
    fped.write("\n")
    if n % 100 == 0:
        print(f"      {n} samples", flush=True)

fped.close()
fmeta.close()
fgeno.close()
print(f"      wrote {n} samples", flush=True)
if bad_tokens:
    print("      WARNING unexpected genotype tokens:", bad_tokens, flush=True)

## 4. duplicate, twin and control pairs from AF31
print("[4/4] copying duplicate/twin/control pairs ...", flush=True)
with open(AF31, newline="") as fh, open(OUT_DUP, "w") as fo:
    r = csv.reader(fh)
    hdr = next(r)
    fo.write("pair\tfirst_sample\tsecond_sample\test_PI\tprocess_control\n")
    for row in r:
        if len(row) < 12 or not row[0].strip():
            continue
        fo.write(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[7]}\t{row[11]}\n")

print("DONE. n_samples=%d n_markers=%d" % (n, n_markers))
