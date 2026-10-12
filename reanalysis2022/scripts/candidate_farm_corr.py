#!/usr/bin/env python3
## correlation of candidate SNP allele frequency and CWD prevalence across farms (>= 10 deer), with p-values,
## BH FDR across candidates and the leave-one-farm-out range
## writes data/candidate_farm_corr.csv and figures/fig2022_04_candidate_structure.png (in reanalysis2022)
import csv, os
import numpy as np
from collections import defaultdict
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

## paths are relative to the reanalysis2022 folder
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(BASE)
MIN_FARM_N = 10

## candidate genes, farms and genotypes
gene = {}
for r in csv.reader(open("data/all_candidates.csv")):
    if len(r) >= 2 and r[1] != "marker":
        gene[r[1].strip()] = r[0].strip()
farm = {r["CMSID"]: r["farm"] for r in csv.DictReader(open("data/cov_1170.csv"))}

rows = list(csv.DictReader(open("data/candidate_genos_1170.tsv"), delimiter="\t"))
markers = [c for c in rows[0] if c.startswith("Affx") or c.startswith("G96S")]
print(f"{len(rows)} deer x {len(markers)} candidate markers")

## CWD prevalence per farm (farms with >= 10 phenotyped deer)
pheno = defaultdict(list)
for r in rows:
    f = farm.get(r["CMSID"], "")
    if f and r["Binary"] in ("0", "1"):
        pheno[f].append(int(r["Binary"]))
farms = {f: (np.mean(v), len(v)) for f, v in pheno.items() if len(v) >= MIN_FARM_N}
print(f"{len(farms)} farms with >= {MIN_FARM_N} phenotyped deer")

## correlation weighted by the number of deer per farm
def wcorr(x, y, w):
    x, y, w = map(np.asarray, (x, y, w)); w = w / w.sum()
    mx, my = (w * x).sum(), (w * y).sum()
    d = np.sqrt((w * (x - mx) ** 2).sum() * (w * (y - my) ** 2).sum())
    return np.nan if d == 0 else (w * (x - mx) * (y - my)).sum() / d

## per candidate (>= 4 farms): r, p, weighted r and leave-one-farm-out range
out = []
for m in markers:
    g = defaultdict(list)
    for r in rows:
        f = farm.get(r["CMSID"], "")
        if f in farms and r[m] in ("0", "1", "2"):
            g[f].append(int(r[m]))
    d = [(f, np.mean(v) / 2, farms[f][0], len(v)) for f, v in g.items() if len(v) >= MIN_FARM_N]
    if len(d) < 4:
        continue
    fr = np.array([t[1] for t in d]); pv = np.array([t[2] for t in d]); n = np.array([t[3] for t in d])
    if fr.std() == 0:
        continue
    r_u, p = stats.pearsonr(fr, pv)
    loo = [np.corrcoef(np.delete(fr, i), np.delete(pv, i))[0, 1] for i in range(len(d))]
    imax = int(np.argmax([abs(l - r_u) for l in loo]))
    out.append(dict(marker=m, gene=gene.get(m, "?"), k_farms=len(d),
                    r=round(r_u, 3), p=round(p, 4), r_nweighted=round(wcorr(fr, pv, n), 3),
                    loo_min=round(min(loo), 3), loo_max=round(max(loo), 3),
                    most_influential_farm=d[imax][0], r_without=round(loo[imax], 3)))

## Benjamini-Hochberg FDR across all candidates tested
ps = np.array([o["p"] for o in out]); order = np.argsort(ps); M = len(ps)
q = np.empty(M); prev = 1.0
for rank, i in enumerate(order[::-1]):
    prev = min(prev, ps[i] * M / (M - rank)); q[i] = prev
for o, qq in zip(out, q):
    o["q_BH"] = round(float(qq), 4)

## write table, sorted by |r|
out.sort(key=lambda o: -abs(o["r"]))
with open("data/candidate_farm_corr.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)

## summary
sig = sum(1 for o in out if o["p"] < 0.05)
sigq = sum(1 for o in out if o["q_BH"] < 0.05)
robust = sum(1 for o in out if o["p"] < 0.05 and np.sign(o["loo_min"]) == np.sign(o["loo_max"])
             and min(abs(o["loo_min"]), abs(o["loo_max"])) > 0.3)
print(f"\ncandidates analysed: {len(out)}")
print(f"  nominally significant (p<0.05) : {sig}/{len(out)} ({sig/len(out)*100:.0f}%)")
print(f"  survive BH-FDR    (q<0.05)     : {sigq}/{len(out)} ({sigq/len(out)*100:.0f}%)")
print(f"  significant AND LOO-robust     : {robust}/{len(out)} ({robust/len(out)*100:.0f}%)")

print(f"\n{'gene':16s}{'k':>4}{'r':>8}{'p':>9}{'q(BH)':>9}{'LOO range':>18}  drop-worst")
for o in out[:14]:
    print(f"{o['gene'][:15]:16s}{o['k_farms']:>4}{o['r']:>8.2f}{o['p']:>9.3f}{o['q_BH']:>9.3f}"
          f"   [{o['loo_min']:+.2f},{o['loo_max']:+.2f}]   {o['most_influential_farm']}→{o['r_without']:+.2f}")

## plot the top 12 candidates by |r|
top = out[:12]
fig, axes = plt.subplots(3, 4, figsize=(17, 10.5))
for ax, o in zip(axes.ravel(), top):
    m = o["marker"]; g = defaultdict(list)
    for r in rows:
        f = farm.get(r["CMSID"], "")
        if f in farms and r[m] in ("0", "1", "2"):
            g[f].append(int(r[m]))
    d = [(f, np.mean(v) / 2, farms[f][0] * 100, len(v)) for f, v in g.items() if len(v) >= MIN_FARM_N]
    fr = np.array([t[1] for t in d]); pv = np.array([t[2] for t in d]); n = np.array([t[3] for t in d])
    ok = o["p"] < 0.05
    ax.scatter(pv, fr, s=np.sqrt(n) * 14, c=["#b2182b" if ok else "#888888"], edgecolor="k", alpha=.85, zorder=3)
    z = np.polyfit(pv, fr, 1); xs = np.linspace(pv.min(), pv.max(), 40)
    ax.plot(xs, np.polyval(z, xs), "--", color="#b2182b" if ok else "#888888", lw=1.4)
    for f, frq, p_, nn in d:
        if f == o["most_influential_farm"]:                # ring the most influential farm
            ax.scatter([p_], [frq], s=np.sqrt(nn) * 14 + 220, facecolors="none",
                       edgecolors="k", linestyle=":", lw=1.4, zorder=4)
    flag = "" if ok else "  n.s."
    ax.set_title(f"{o['gene'][:22]}\nr={o['r']:+.2f}  p={o['p']:.3f}  q={o['q_BH']:.2f}{flag}",
                 fontsize=9.5, color="#b2182b" if ok else "#555555",
                 fontweight="bold" if ok else "normal")
    ax.text(0.03, 0.96, f"drop {o['most_influential_farm']}→r={o['r_without']:+.2f}",
            transform=ax.transAxes, va="top", fontsize=7.2, style="italic",
            bbox=dict(boxstyle="round,pad=0.25", fc="#fff6d5", ec="#d8c88a"))
    ax.tick_params(labelsize=7.5)
for ax in axes.ravel()[len(top):]:
    ax.axis("off")
fig.suptitle("Candidate loci: farm allele-frequency vs farm CWD-prevalence (top 12 of "
             f"{len(out)} candidates by |r|)\n"
             f"only {sig}/{len(out)} candidates reach p<0.05 and {sigq}/{len(out)} survive FDR — "
             "each correlation rests on ~16 farms and is dominated by one herd (dotted ring)",
             fontsize=13)
fig.supxlabel("farm CWD prevalence (%)", fontsize=11)
fig.supylabel("candidate minor-allele frequency", fontsize=11)
fig.tight_layout(rect=[0.015, 0.015, 1, 0.93])
fig.savefig("figures/fig2022_04_candidate_structure.png", dpi=150)
print("\nwrote data/candidate_farm_corr.csv and figures/fig2022_04_candidate_structure.png")
