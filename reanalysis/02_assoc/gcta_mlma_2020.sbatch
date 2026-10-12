#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --partition=acpu
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --account=YOUR_ACCOUNT          # set to your allocation
#SBATCH --job-name=gcta_mlma_loco
#SBATCH --error=logs/gcta_mlma_loco-%A-%a.err.txt 
#SBATCH --output=logs/gcta_mlma_loco-%A-%a.out.txt

## GCTA --mlma and --mlma-loco for the 2020 dataset (Table 1, 2020 GCTA rows)
## deer_bin, deer_scr = reanalysis/01_structure/deer_auto (chr 1-29, MAF >= 0.05, 120,113 SNPs) with the binary or tissue score phenotype
## pheno_*.gcta = FID IID phenotype
source ~/.bashrc

## Activate an environment with gcta64 (GCTA 1.94.1)

## Set your master path (folder with deer_bin, deer_scr and pheno_*.gcta)
MASTER=<dir>/Seabury_data/reanalysis/02_assoc

cd $MASTER

## mlma run
gcta64 --bfile deer_bin --autosome --autosome-num 29 --make-grm --out gcta_grm --thread-num 4
gcta64 --mlma --bfile deer_bin --grm gcta_grm --pheno pheno_bin.gcta --autosome-num 29 --out gcta_mlma_bin --thread-num 4
gcta64 --mlma --bfile deer_scr --grm gcta_grm --pheno pheno_scr.gcta --autosome-num 29 --out gcta_mlma_scr --thread-num 4

## loco
gcta64 --mlma-loco --bfile deer_bin --grm gcta_grm --pheno pheno_bin.gcta --autosome-num 29 --out gcta_loco_bin --thread-num 4
gcta64 --mlma-loco --bfile deer_scr --grm gcta_grm --pheno pheno_scr.gcta --autosome-num 29 --out gcta__loco_scr --thread-num 4
