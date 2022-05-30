#!/bin/sh
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=80:00:00
#SBATCH --job-name=fit_sample_type
#SBATCH -o launch_scripts/cluster_logs/fit_sample_type_%a_%A.log 
#SBATCH --mem-per-cpu=7900MB
#SBATCH --array=1-24

# -- commands you want to execute -- #

patient=`head -n ${SLURM_ARRAY_TASK_ID} analysis/ATAC/CLAIRE/Epigenetic/input_data/patients.txt | tail -1`

Rscript analysis/ATAC/CLAIRE/Epigenetic/04-ATAC_regional_recurrent_peaks_deseq_LRT.R \
${patient} ~sample_type+purity ~sample_type effect_of_purity_only_nfreads_LRT
