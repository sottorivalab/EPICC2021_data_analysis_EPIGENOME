#!/bin/sh
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=00:10:00
#SBATCH --job-name=add_lp_samples
#SBATCH -o launch_scripts/cluster_logs/add_lp_samples_%a_%j.log
#SBATCH --array=1-240
#SBATCH --mem-per-cpu=1000MB

set -eo pipefail
module load R/4.0.5
TMPDIR=$(mktemp -d /tmp/$USER__$SLURM_JOBID-XXXX)

Rscript  analysis/WGS/03-TREES/2.1-assign_lowpass_samples_to_tree.R "$SLURM_ARRAY_TASK_ID"

echo "Tmp size is (after job run): $(du -sh $tempdir)"
rm -rf $tempdir

