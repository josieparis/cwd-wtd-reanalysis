#!/bin/bash

## PLINK steps for the 2022 dataset: GRM, GRM without chr13 (PRNP) and PCA

## Set your master path
MASTER=<dir>/Seabury_data
cd $MASTER

## make the plink files from the Dryad genotypes first (120,020 autosomal SNPs):
# unzip -p doi_10_5061_dryad_wh70rxwnt__v20220622.zip 1170_WTD_GWAAs_IncludingGxE.tsv | python3 reanalysis2022/scripts/build_bed_2022.py reanalysis2022/work/deer2022
# plink --tfile reanalysis2022/work/deer2022 --chr-set 29 --make-bed --out reanalysis2022/work/deer2022

WORK_DIR=reanalysis2022/work
B=$WORK_DIR/deer2022
S=reanalysis2022/01_structure

## GRM (VanRaden), autosomes, geno 0.05, maf 0.01 (119,962 SNPs)
plink --bfile $B --chr 1-29 --chr-set 29 --geno 0.05 --maf 0.01 --make-rel square gz --out $S/grm2022

## same GRM without chr13 (115,272 SNPs), for the polygenic remainder in Fig 3
plink --bfile $B --chr 1-29 --chr-set 29 --geno 0.05 --maf 0.01 --not-chr 13 --make-rel square gz --out $S/grm2022_no13

## PCA (Fig 1a)
plink --bfile $B --chr 1-29 --chr-set 29 --geno 0.05 --maf 0.01 --pca 20 --out $S/pca2022
