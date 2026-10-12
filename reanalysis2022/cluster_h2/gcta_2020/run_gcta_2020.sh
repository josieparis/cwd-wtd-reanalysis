#!/bin/bash

## GCTA-GREML heritability for the 2020 dataset (n = 807), with and without covariates
## run inside the gcta_2020 folder: GCTA=$(which gcta64) bash run_gcta_2020.sh

set -e
GCTA=${GCTA:-gcta64}

## binary
$GCTA --grm deer2020 --pheno pheno_bin.txt                                                       --reml --out b_none
$GCTA --grm deer2020 --pheno pheno_bin.txt --covar covar_region.txt --qcovar qcovar_age.txt      --reml --out b_region
$GCTA --grm deer2020 --pheno pheno_bin.txt --covar covar_region.txt --qcovar qcovar_agescore.txt --reml --out b_region_score
## score (CWD tissue score)
$GCTA --grm deer2020 --pheno pheno_scr.txt                                                        --reml --out s_none
$GCTA --grm deer2020 --pheno pheno_scr.txt --covar covar_region.txt --qcovar qcovar_age.txt       --reml --out s_region

## print h2
echo; echo "============ 2020 GCTA h2 (V(G)/Vp) ============"
for f in b_none b_region b_region_score s_none s_region; do
  printf "%-16s " "$f"; grep "V(G)/Vp" $f.hsq | head -1
done
