#!/bin/bash

## GEMMA SNP heritability (pve) for the 2022 dataset
## GEMMA builds its own centred GRM, so this is an independent check of the GCTA estimates
## run inside the gemma folder: GEMMA=/path/to/gemma bash run_gemma_h2.sh (or module load gemma)

## inputs
## geno.bimbam.gz: 120,020 autosomal SNPs x 1171 deer, BIMBAM mean genotypes
## pheno.txt: binary CWD (1170 phenotyped), pheno_region.txt: complete cases (n = 1166),
## pheno_farm.txt: complete cases on farms with >= 10 deer (n = 1151), pheno_grm.txt: all 1171 (dummy, for the GRM)
## covar_intercept.txt / covar_region.txt / covar_farm.txt (first column = intercept)
## the GRM is built on all 1171 deer so it matches every -lmm run, and each run drops its own missing phenotypes

set -e
GEMMA=${GEMMA:-gemma}
gunzip -kf geno.bimbam.gz

## GEMMA centred relatedness matrix (all 1171 deer)
$GEMMA -g geno.bimbam -p pheno_grm.txt -gk 1 -o grm2022        # writes output/grm2022.cXX.txt

K=output/grm2022.cXX.txt
## pve estimate in the null model = SNP heritability, for each covariate set
$GEMMA -g geno.bimbam -p pheno.txt        -k $K -c covar_intercept.txt -lmm 1 -n 1 -o pve_nocov
$GEMMA -g geno.bimbam -p pheno_region.txt -k $K -c covar_region.txt    -lmm 1 -n 1 -o pve_region
$GEMMA -g geno.bimbam -p pheno_farm.txt   -k $K -c covar_farm.txt      -lmm 1 -n 1 -o pve_farm

## print results
echo; echo "================ RESULTS ================"
grep -H "pve estimate in the null model" output/pve_*.log.txt
grep -H "se(pve)"                          output/pve_*.log.txt
