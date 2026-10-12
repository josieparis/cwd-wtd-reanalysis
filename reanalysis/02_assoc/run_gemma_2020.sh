#!/bin/bash

## GEMMA 0.98.5 single-GRM scans for the 2020 dataset, as run (commands from the GEMMA logs in reanalysis/02_assoc/output/*.log.txt)
## run inside reanalysis/02_assoc

## inputs
## deer_bin, deer_scr = reanalysis/01_structure/deer_auto (chr 1-29, MAF >= 0.05, 120,113 SNPs) with binary CWD status
## or tissue score (0/1/2) as the .fam phenotype (from 00_data/samples.tsv)
## deer_wX_bin, deer_wX_scr = the same plus chr X (123,236 SNPs), to match the authors' Table S1 which lists X-linked markers:
# plink --bfile ../00_data/deer --chr 1-30 --chr-set 29 --geno 0.05 --maf 0.05 --make-bed --out deer_wX
## gemma_grm.txt = reanalysis/01_structure/grm.rel.gz (plink --make-rel, autosomes) unzipped, used for every scan
## covar_base.txt = intercept, sex, age, Midwest, South (Northeast is the reference)
## covar_score.txt = covar_base.txt + tissue score (the authors' "+ CWD score" model)

## gemma executable
GEMMA=gemma

## Table 1 (2020 GEMMA rows) and threshold counts (13_threshold_counts.py)
$GEMMA -bfile deer_bin -k gemma_grm.txt -lmm 4 -o gemma_bin
$GEMMA -bfile deer_scr -k gemma_grm.txt -lmm 4 -o gemma_scr

## the authors' five Table S1 models (Table S1 recount, tissue score analysis, Supplementary Fig 2d)
$GEMMA -bfile deer_wX_bin -k gemma_grm.txt                     -lmm 4 -o gwaa_bin_nocov
$GEMMA -bfile deer_wX_bin -k gemma_grm.txt -c covar_base.txt   -lmm 4 -o gwaa_bin_cov
$GEMMA -bfile deer_wX_bin -k gemma_grm.txt -c covar_score.txt  -lmm 4 -o gwaa_bin_covscore
$GEMMA -bfile deer_wX_scr -k gemma_grm.txt                     -lmm 4 -o gwaa_scr_nocov
$GEMMA -bfile deer_wX_scr -k gemma_grm.txt -c covar_base.txt   -lmm 4 -o gwaa_scr_cov

## GEMMA writes to ./output/
## gemma_{bin,scr}_results.tsv (read by 13_threshold_counts.py) are the rs, beta, se, p_wald, p_lrt and p_score
## columns of output/gemma_{bin,scr}.assoc.txt, with rs renamed to marker
