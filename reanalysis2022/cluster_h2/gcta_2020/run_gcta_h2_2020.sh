#!/bin/bash
#SBATCH --job-name=gcta_h2_2020
#SBATCH --time=00:30:00
#SBATCH --mem=15360M
#SBATCH --cpus-per-task=4
#SBATCH --account=YOUR_ACCOUNT          # set to your allocation
#SBATCH --partition=acpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --qos=normal
#SBATCH --output=logs/gcta_h2_2020_%A.out
#SBATCH --error=logs/gcta_h2_2020_%A.err

## GCTA-GREML heritability for the 2020 dataset (n = 807), as in the 2022 run (run_gcta_h2.sh)
## no farm IDs in 2020, so the models are no covariates, + sex/age/region and + sex/age/region + tissue score
## GRM = deer2020.grm.* (our VanRaden GRM), cases 284/807 = 0.352; submit from inside the gcta_2020 folder
set -e
GCTA=${GCTA:-gcta64}                    # GCTA 1.94.1
$GCTA --version 2>&1 | head -3

## observed scale
$GCTA --grm deer2020 --pheno pheno_bin.txt                                                       --reml --out b_none
$GCTA --grm deer2020 --pheno pheno_bin.txt --covar covar_region.txt --qcovar qcovar_age.txt      --reml --out b_region
$GCTA --grm deer2020 --pheno pheno_bin.txt --covar covar_region.txt --qcovar qcovar_agescore.txt --reml --out b_region_score
$GCTA --grm deer2020 --pheno pheno_scr.txt                                                        --reml --out s_none
$GCTA --grm deer2020 --pheno pheno_scr.txt --covar covar_region.txt --qcovar qcovar_age.txt       --reml --out s_region

## liability scale (binary), GCTA --prevalence as in Seabury et al.
for K in 0.05 0.10 0.26; do
  $GCTA --grm deer2020 --pheno pheno_bin.txt --reml --prevalence $K --out b_none_liab_$K
done

## print results
echo; echo "================ 2020 GCTA h2 (V(G)/Vp) ================"
for f in b_none b_region b_region_score s_none s_region; do
  echo "--- $f (observed) ---"; grep -E "V\(G\)/Vp|Pval|^n " $f.hsq 2>/dev/null
done
echo "--- liability scale, no covariates ---"
for K in 0.05 0.10 0.26; do
  echo "prevalence=$K:"; grep -E "V\(G\)/Vp_L" b_none_liab_$K.hsq 2>/dev/null
done
