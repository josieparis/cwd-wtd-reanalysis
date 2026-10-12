#!/bin/bash
#SBATCH --job-name=gcta_mlma_2022_final
#SBATCH --time=00:45:00
#SBATCH --mem=15360M
#SBATCH --cpus-per-task=4
#SBATCH --account=YOUR_ACCOUNT          # set to your allocation
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=logs/gcta_mlma_2022_final_%A.out

## final 2022 association scans for Table 1, submit from inside the mlma_results folder
## --mlma on the PLINK --make-rel GRM, and --mlma-loco with --autosome-num 29 (the default of 22 drops chr 23-29)
set -e
cd $SLURM_SUBMIT_DIR
GCTA=${GCTA:-gcta64}                    # GCTA 1.94.1
GRM=../cluster_h2/gcta_plink/deer2022p

## run association scans
$GCTA --bfile deer2022 --grm $GRM --pheno pheno_bin.txt --autosome-num 29 --thread-num 4 --mlma      --out mlma_bin_plink
$GCTA --bfile deer2022            --pheno pheno_bin.txt --autosome-num 29 --thread-num 4 --mlma-loco --out loco_bin_chr29

## number of SNPs and chromosomes, hits at P < 5e-8 and P < 5e-5, and the PRNP line
for F in mlma_bin_plink.mlma loco_bin_chr29.loco.mlma; do
  echo "--- $F: $(($(wc -l < $F)-1)) SNPs on $(awk "NR>1{print \$1}" $F | sort -u | wc -l) chromosomes"
  awk "NR>1{p=\$NF; if(p<5e-8)c8++; if(p<5e-5)c5++} END{printf \"P<5e-8=%d  P<5e-5=%d\n\",c8,c5}" $F
  awk "NR==1 || \$2==\"G96S_G_to_A_at286\"" $F
done
