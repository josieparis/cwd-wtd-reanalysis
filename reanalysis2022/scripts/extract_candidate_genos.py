#!/usr/bin/env python3
## keep only the candidate marker columns (reanalysis2022/data/all_candidates.csv) of the 2022 genotype file on stdin
## writes reanalysis2022/data/candidate_genos_1170.tsv and candidate_genos_missing.txt (candidates not genotyped)
## unzip -p doi_10_5061_dryad_wh70rxwnt__v20220622.zip 1170_WTD_GWAAs_IncludingGxE.tsv | python3 reanalysis2022/scripts/extract_candidate_genos.py
import sys, csv, os

## run from the main project folder
BASE = "reanalysis2022"

## candidate markers and genes
cands = {}
with open(os.path.join(BASE, "data/all_candidates.csv")) as fh:
    for r in csv.DictReader(fh):
        if r.get("marker"):
            cands[r["marker"].strip()] = r["gene"].strip()

## candidate columns in the genotype file
hdr = sys.stdin.readline().rstrip("\n").split("\t")
meta = hdr[:4]                                    # CMSID, GeneralRegion, Sex, Binary
idx = [i for i, m in enumerate(hdr) if i >= 4 and m in cands]
keep = [hdr[i] for i in idx]
missing = sorted(set(cands) - set(keep))
print(f"[extract] {len(cands)} candidates requested; {len(keep)} found in genotypes; "
      f"{len(missing)} absent", file=sys.stderr, flush=True)

## write candidate genotypes
out_path = os.path.join(BASE, "data/candidate_genos_1170.tsv")
n = 0
with open(out_path, "w") as out:
    out.write("\t".join(meta + keep) + "\n")
    for line in sys.stdin:
        f = line.rstrip("\n").split("\t")
        out.write("\t".join(f[:4] + [f[i] for i in idx]) + "\n")
        n += 1
        if n % 200 == 0:
            print(f"[extract] {n} deer", file=sys.stderr, flush=True)

## write candidates missing from the genotypes
with open(os.path.join(BASE, "data/candidate_genos_missing.txt"), "w") as mf:
    mf.write("\n".join(f"{m}\t{cands[m]}" for m in missing) + "\n")
print(f"[extract] wrote {out_path}  ({n} deer x {len(keep)} candidate markers)", file=sys.stderr)
