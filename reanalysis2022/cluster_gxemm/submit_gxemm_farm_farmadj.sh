#!/bin/bash
#SBATCH --job-name=gxemm_farm_redo
#SBATCH --time=08:00:00
#SBATCH --mem=15360M
#SBATCH --cpus-per-task=4
#SBATCH --account=YOUR_ACCOUNT          # set to your allocation
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1

## re-run only the farm and farm-adjusted GxEMM models (the region run had already finished on the PLINK GRM)
## submit from inside the cluster_gxemm folder
cd "$SLURM_SUBMIT_DIR"
RS=${RS:-Rscript}                       # R 4.5.3 with GxEMM 1.0 and msm installed
chmod +x ldak5.linux

## run models
echo "===== FARM     $(date) ====="; "$RS" run_gxemm_farm_cluster.R    || echo ">>> farm FAILED"
echo "===== FARMADJ  $(date) ====="; "$RS" run_gxemm_farmadj_cluster.R || echo ">>> farmadj FAILED"
echo "===== DONE FARM+FARMADJ $(date) ====="
