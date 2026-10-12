#!/bin/bash

## PLINK 1.9 steps for the 2020 dataset (n = 807), as run (see reanalysis/01_structure/*.log)
## run after 00_build_dataset.py has written reanalysis/00_data/deer.map and reanalysis/work/deer.ped

## Set your master path
MASTER=<dir>/Seabury_data
cd $MASTER

## large intermediates go in WORK_DIR, kept out of the repository
WORK_DIR=reanalysis/work; mkdir -p "$WORK_DIR"
D=reanalysis/00_data
S=reanalysis/01_structure

## 0) binary plink files from the Dryad genotypes (00_build_dataset.py)
plink --ped "$WORK_DIR/deer.ped" --map $D/deer.map --chr-set 29 --make-bed --out $D/deer

## 1) autosomal association set: chr 1-29, geno 0.05, maf 0.05 (120,113 SNPs)
plink --bfile $D/deer --chr 1-29 --chr-set 29 --geno 0.05 --maf 0.05 --make-bed --out $S/deer_auto

## 2) LD pruning (--indep-pairwise 200 25 0.2) and 3) PCA on the pruned set (Fig 1a)
plink --bfile $S/deer_auto --chr-set 29 --indep-pairwise 200 25 0.2 --out $S/prune
plink --bfile $S/deer_auto --chr-set 29 --extract $S/prune.prune.in --pca 20 header tabs --out $S/pca

## 4) GRM (VanRaden, --make-rel), autosomes, geno 0.05, maf 0.01, not pruned (120,368 SNPs)
plink --bfile $D/deer --chr 1-29 --chr-set 29 --geno 0.05 --maf 0.01 --make-rel square gz --out $S/grm
## same GRM without chr13 (PRNP), for the polygenic remainder and variance partition models
plink --bfile $D/deer --chr 1-29 --chr-set 29 --geno 0.05 --maf 0.01 --not-chr 13 --make-rel square gz --out $S/grm_no13

## 5) IBD (PI_HAT) on the pruned set, for relative pair counts and the duplicate check
plink --bfile $S/deer_auto --chr-set 29 --extract $S/prune.prune.in --genome --out $S/ibd

## 6) PRNP marker QC and additive dosages (--recode A) for the association and prediction scripts
## the two .raw files were not archived, these commands remake them
## deer_auto.raw is read by 04_gwas_lambda.py and 07_mlma_loco.py, prnp.raw by 06, 09 and 11
## (G96S_G_to_A_at286_A = dosage of the serine allele)
plink --bfile $D/deer --chr-set 29 --snps G96S_G_to_A_at286,Q226K_C_to_A_at676,NewMutationCMS1 --freq --missing --out $S/prnp_qc
plink --bfile $D/deer --chr-set 29 --snp G96S_G_to_A_at286 --recode A --out "$WORK_DIR/prnp"       # prnp.raw
plink --bfile $S/deer_auto --chr-set 29 --recode A --out "$WORK_DIR/deer_auto"                     # deer_auto.raw
