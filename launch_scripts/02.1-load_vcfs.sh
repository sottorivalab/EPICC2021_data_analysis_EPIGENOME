#!/bin/sh
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=02:30:00
#SBATCH --job-name=ls_load_vcf
#SBATCH -o launch_scripts/cluster_logs/load_vcfs_%a_%j.log
#SBATCH --array=1-60
#SBATCH --mem-per-cpu=2500MB

set -eo pipefail
module load R/4.0.5
TMPDIR=$(mktemp -d /tmp/$USER__$SLURM_JOBID-XXXX)

#Start the job here.
Rscript analysis/WGS/02-SNV/1-load_snv_data.R "$SLURM_ARRAY_TASK_ID"

echo "Tmp size is (after job run): $(du -sh $tempdir)"
rm -rf $tempdir

