#!/bin/bash
#SBATCH --job-name=gxemm_cwd_plink
#SBATCH --time=12:00:00
#SBATCH --mem=15360M
#SBATCH --cpus-per-task=4
#SBATCH --account=YOUR_ACCOUNT          # set to your allocation
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1

## run the three 2022 GxEMM models (region, farm, farm-adjusted) on the PLINK GRM in one job
## submit from inside the cluster_gxemm folder; uses the R in the gatk-tools conda env
## mem = 4 x 3840 MB, so acpu bills exactly 4 cores
set -o pipefail
cd "$SLURM_SUBMIT_DIR"
RS=${RS:-Rscript}                       # R 4.5.3 with GxEMM 1.0 and msm installed
chmod +x ldak5.linux

## run models
echo "===== REGION   $(date) ====="; "$RS" run_gxemm_cluster.R      || echo ">>> region FAILED"
echo "===== FARM     $(date) ====="; "$RS" run_gxemm_farm_cluster.R   || echo ">>> farm FAILED"
echo "===== FARMADJ  $(date) ====="; "$RS" run_gxemm_farmadj_cluster.R || echo ">>> farmadj FAILED"
echo "===== ALL DONE $(date) ====="
