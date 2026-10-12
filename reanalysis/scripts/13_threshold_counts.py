#!/usr/bin/env python3
## count markers at the authors' polygenic threshold (P <= 5e-5), and at 5e-7 and 5e-8, for each model
## also counts independent loci (SNPs more than 1 Mb apart) so one signal isn't counted many times
## writes reanalysis/01_structure/threshold_counts.csv
import gzip, csv, os

PRNP = "G96S_G_to_A_at286"
## paths are relative to the reanalysis folder
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(BASE)

## SNP positions
pos = {}
for line in open("02_assoc/deer_bin.bim"):
    f = line.split()
    pos[f[1]] = (f[0], int(f[3]))

def read_tsv(path, key, pcol, gz=False):
    op = gzip.open(path, "rt") if gz else open(path)
    with op as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        return {r[key]: float(r[pcol]) for r in rd
                if r.get(pcol) not in (None, "", "NA", "nan")}

## p-values for each model
M = {}
for trait, suf in [("binary", "bin"), ("score", "scr")]:
    loco = f"02_assoc/output/gcta_loco_{suf}.loco.mlma"
    if not os.path.exists(loco):
        loco = f"02_assoc/output/gcta__loco_{suf}.loco.mlma"
    M[trait] = {
        "Naive OLS (no GRM)":        read_tsv(f"01_structure/gwas_{trait}.tsv.gz", "marker", "p_naive", gz=True),
        "OLS + 10 genomic PCs":      read_tsv(f"01_structure/gwas_{trait}.tsv.gz", "marker", "p_pc",    gz=True),
        "GEMMA (single GRM, LMM)":   read_tsv(f"02_assoc/gemma_{suf}_results.tsv", "marker", "p_wald"),
        "GCTA --mlma (single GRM)":  read_tsv(f"02_assoc/output/gcta_mlma_{suf}.mlma", "SNP", "p"),
        "GCTA --mlma-loco":          read_tsv(loco, "SNP", "p"),
    }

common = set.intersection(*[set(M["binary"][k]) for k in M["binary"]])
print(f"common markers tested across all models: {len(common):,}")

def indep_loci(hits, window=1_000_000):
    # keep the best SNP, drop others within 1 Mb on the same chr, repeat
    taken = []
    for m in sorted(hits, key=lambda x: hits[x]):
        if m not in pos:
            taken.append(m); continue
        c, b = pos[m]
        if all(not (pos[t][0] == c and abs(pos[t][1] - b) < window) for t in taken if t in pos):
            taken.append(m)
    return len(taken)

## counts per model
rows = []
for trait in ("binary", "score"):
    for name, d in M[trait].items():
        dd = {m: v for m, v in d.items() if m in common}
        hits5 = {m: v for m, v in dd.items() if v <= 5e-5}
        rows.append(dict(
            trait=trait, model=name,
            n_5e05=len(hits5), loci_5e05=indep_loci(hits5),
            n_5e07=sum(1 for v in dd.values() if v <= 5e-7),
            n_5e08=sum(1 for v in dd.values() if v <= 5e-8),
            PRNP_p=f"{dd.get(PRNP):.3e}" if dd.get(PRNP) else "NA",
        ))

## write table
with open("01_structure/threshold_counts.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)

for trait in ("binary", "score"):
    print(f"\n=== {trait} ===")
    print(f"{'model':30s}{'@5e-5':>9s}{'loci':>7s}{'@5e-7':>9s}{'@5e-8':>9s}{'PRNP p':>12s}")
    for r in rows:
        if r["trait"] == trait:
            print(f"{r['model']:30s}{r['n_5e05']:>9,}{r['loci_5e05']:>7,}{r['n_5e07']:>9,}{r['n_5e08']:>9,}{r['PRNP_p']:>12s}")
print("\nwrote 01_structure/threshold_counts.csv")
