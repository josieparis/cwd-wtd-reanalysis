#!/bin/bash
#SBATCH --job-name=gcta_h2_plink
#SBATCH --time=00:30:00
#SBATCH --mem=15360M
#SBATCH --cpus-per-task=4
#SBATCH --account=YOUR_ACCOUNT          # set to your allocation
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=gcta_h2_plink_%j.out

## GCTA-GREML heritability for the 2022 dataset on the PLINK GRM (deer2022p.grm.*), observed and liability scale
## each model uses the original paper's complete-case sample (1,170 / 1,166 / 1,151)
## inputs from build_gcta_plink_inputs.R; submit from inside the gcta_plink folder: sbatch run_gcta_plink.sh
set -e
cd "${SLURM_SUBMIT_DIR:-.}"
GCTA=${GCTA:-gcta64}                    # GCTA 1.94.1
T="--thread-num ${SLURM_CPUS_PER_TASK:-4}"
mkdir -p results

## observed scale
$GCTA --grm deer2022p --pheno pheno.txt                                                    --reml $T --out results/h2_nocov
$GCTA --grm deer2022p --pheno pheno.txt --covar covar_region.txt --qcovar qcovar_age.txt   --reml $T --out results/h2_region
$GCTA --grm deer2022p --pheno pheno.txt --keep keep_farm10.txt \
      --covar covar_farm10.txt --qcovar qcovar_age.txt                                     --reml $T --out results/h2_farm

## liability scale, no covariates (GCTA --prevalence, Lee et al. 2011, as in Seabury et al.)
for K in 0.05 0.10 0.26; do
  $GCTA --grm deer2022p --pheno pheno.txt --reml --prevalence $K $T --out results/h2_nocov_liab_$K
done

## print results
echo; echo "================ RESULTS (PLINK GRM) ================"
for f in h2_nocov h2_region h2_farm; do
  printf "%-10s " $f; awk -F'\t' '$1=="V(G)/Vp"{printf "h2 = %s +- %s  ",$2,$3} $1=="n"{printf "n = %s",$2}' results/$f.hsq; echo
done
for K in 0.05 0.10 0.26; do
  printf "liab K=%s  " $K; awk -F'\t' '$1=="V(G)/Vp_L"{printf "h2L = %s +- %s",$2,$3}' results/h2_nocov_liab_$K.hsq; echo
done
