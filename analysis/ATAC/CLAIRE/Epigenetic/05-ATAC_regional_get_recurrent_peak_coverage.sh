#!/bin/sh
#SBATCH -n 6
#SBATCH -N 1
#SBATCH --time=80:00:00
#SBATCH --job-name=fit_region_purity
#SBATCH -o launch_scripts/cluster_logs/get_peak_covrage_%a_%A.log 
#SBATCH --mem-per-cpu=1500MB
#SBATCH --array=1-160

Rscript analysis/ATAC/CLAIRE/Epigenetic/05-ATAC_regional_get_recurrent_peak_coverage.R ${SLURM_ARRAY_TASK_ID}
