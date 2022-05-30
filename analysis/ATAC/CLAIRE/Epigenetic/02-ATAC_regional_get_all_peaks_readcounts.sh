#!/bin/bash
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --time=80:00:00
#SBATCH --job-name=megabulks
#SBATCH -o launch_scripts/cluster_logs/get_readcounts_%a.log 
#SBATCH --mem-per-cpu=7900MB
#SBATCH --array=1-24

patient=`head -n ${SLURM_ARRAY_TASK_ID} analysis/ATAC/CLAIRE/Epigenetic/input_data/patients.txt | tail -1`
Rscript analysis/ATAC/CLAIRE/Epigenetic/02-ATAC_regional_get_recurrent_peaks_readcounts.R -p "$patient" -r analysis/ATAC/CLAIRE/Epigenetic/datasets/all_peaks.rds 
