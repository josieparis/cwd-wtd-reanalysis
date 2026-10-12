#!/usr/bin/env python3
## Step 1a: PCA structure, variance explained, PCs vs region and phenotype, and PCA plots
## writes reanalysis/figures/fig1_pca_by_region.png, fig2_pca_by_phenotype.png, fig3_scree.png and 01_structure/pca_merged.tsv
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

## run from the main project folder
ROOT = "reanalysis"
S=f"{ROOT}/01_structure"; F=f"{ROOT}/figures"

## read in PCA and sample info
ev = np.loadtxt(f"{S}/pca.eigenval")
vexp = ev/ev.sum()*100
vec = pd.read_csv(f"{S}/pca.eigenvec", sep="\t")
vec = vec.rename(columns={vec.columns[0]:"FID", vec.columns[1]:"IID"})
meta = pd.read_csv(f"{ROOT}/00_data/samples.tsv", sep="\t", dtype={"clean_id":str})
df = vec.merge(meta, left_on="IID", right_on="clean_id")
pcs=[f"PC{i}" for i in range(1,21)]

## variance explained
print("=== Variance explained by leading PCs ===")
for i in range(10):
    print(f"  PC{i+1}: eigenvalue={ev[i]:8.3f}  var={vexp[i]:6.2f}%  cum={vexp[:i+1].sum():6.2f}%")

## PCs vs region (Kruskal-Wallis) and phenotype (point-biserial, Spearman)
print("\n=== PC ~ REGION (Kruskal-Wallis) & PC ~ phenotype (point-biserial / Spearman) ===")
print(f"{'PC':>4} {'var%':>6} {'KW_region_p':>12} {'r_binary':>9} {'p_binary':>10} {'rho_score':>9} {'p_score':>10}")
reg_groups=[df[df.region==r] for r in ["South","Midwest","East"]]
for i in range(10):
    pc=f"PC{i+1}"
    kw=stats.kruskal(*[g[pc] for g in reg_groups])
    rb,pb=stats.pointbiserialr(df["binary"], df[pc])
    rs,ps=stats.spearmanr(df["score"], df[pc])
    print(f"{pc:>4} {vexp[i]:6.2f} {kw.pvalue:12.2e} {rb:9.3f} {pb:10.2e} {rs:9.3f} {ps:10.2e}")

## region centroids
print("\n=== Region centroids (PC1, PC2) ===")
print(df.groupby("region")[["PC1","PC2","PC3"]].mean().round(4).to_string())

## plot
reg_col={"South":"#1f77b4","Midwest":"#d62728","East":"#2ca02c"}
def scatter(ax,x,y,groups,cmap,title,leg):
    for g,c in cmap.items():
        s=df[df[groups]==g]
        ax.scatter(s[x],s[y],s=14,alpha=.7,c=c,label=f"{g} (n={len(s)})",edgecolors="none")
    ax.set_xlabel(f"{x} ({vexp[int(x[2:])-1]:.1f}%)")
    ax.set_ylabel(f"{y} ({vexp[int(y[2:])-1]:.1f}%)")
    ax.set_title(title);
    if leg: ax.legend(fontsize=7,markerscale=1.3,frameon=False)

## fig 1: PCA by region (PC1-2 and PC3-4)
fig,ax=plt.subplots(1,2,figsize=(12,5))
scatter(ax[0],"PC1","PC2","region",reg_col,"Genomic PCA by region",True)
scatter(ax[1],"PC3","PC4","region",reg_col,"PC3 vs PC4 by region",True)
fig.tight_layout(); fig.savefig(f"{F}/fig1_pca_by_region.png",dpi=150); plt.close(fig)

## fig 2: PCA coloured by phenotype (binary and score)
fig,ax=plt.subplots(1,2,figsize=(12,5))
bcol={0:"#bbbbbb",1:"#e41a1c"}
for b,c in bcol.items():
    s=df[df.binary==b]; ax[0].scatter(s.PC1,s.PC2,s=14,alpha=.7,c=c,
        label=("CWD+ " if b else "non-detect ")+f"(n={len(s)})",edgecolors="none")
ax[0].set_xlabel(f"PC1 ({vexp[0]:.1f}%)");ax[0].set_ylabel(f"PC2 ({vexp[1]:.1f}%)")
ax[0].set_title("PCA by CWD binary status");ax[0].legend(fontsize=7,frameon=False)
scol={0:"#bbbbbb",1:"#fdae61",2:"#d7191c"}
for sc,c in scol.items():
    s=df[df.score==sc]; ax[1].scatter(s.PC1,s.PC2,s=14,alpha=.7,c=c,
        label=f"score {sc} (n={len(s)})",edgecolors="none")
ax[1].set_xlabel(f"PC1 ({vexp[0]:.1f}%)");ax[1].set_ylabel(f"PC2 ({vexp[1]:.1f}%)")
ax[1].set_title("PCA by CWD tissue score");ax[1].legend(fontsize=7,frameon=False)
fig.tight_layout(); fig.savefig(f"{F}/fig2_pca_by_phenotype.png",dpi=150); plt.close(fig)

## fig 3: scree plot
fig,ax=plt.subplots(figsize=(6,4))
ax.bar(range(1,21),vexp,color="#555");ax.set_xlabel("PC");ax.set_ylabel("% variance")
ax.set_title("Scree plot (autosomal, LD-pruned; 99,433 SNPs)")
fig.tight_layout(); fig.savefig(f"{F}/fig3_scree.png",dpi=150); plt.close(fig)

## save merged table
df.to_csv(f"{S}/pca_merged.tsv",sep="\t",index=False)
print("\nfigures written: fig1_pca_by_region, fig2_pca_by_phenotype, fig3_scree")
