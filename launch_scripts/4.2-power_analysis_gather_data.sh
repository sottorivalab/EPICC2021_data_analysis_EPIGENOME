#!/bin/sh
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=00:30:00
#SBATCH --job-name=power_analysis
#SBATCH -o launch_scripts/cluster_logs/power_analysis_%a_%A.log 
#SBATCH --mem-per-cpu=7500MB
#SBATCH --array=1-109

Rscript analysis/ATAC/MEGABULKS/4-power_analysis_gather_data.R ${SLURM_ARRAY_TASK_ID}
