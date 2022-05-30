#!/bin/sh
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=03:00:00
#SBATCH --job-name=lift_for_dnds
#SBATCH -o launch_scripts/cluster_logs/lift_for_dnds_%a_%j.log
#SBATCH --array=1-60
#SBATCH --mem-per-cpu=2500MB

set -eo pipefail
module load R/4.0.5
TMPDIR=$(mktemp -d /tmp/$USER__$SLURM_JOBID-XXXX)

#Start the job here.
Rscript analysis/WGS/02-SNV/5.1-dnds.R "$SLURM_ARRAY_TASK_ID"

echo "Tmp size is (after job run): $(du -sh $tempdir)"
rm -rf $tempdir

