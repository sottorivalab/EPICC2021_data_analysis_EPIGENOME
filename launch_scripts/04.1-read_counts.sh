#!/bin/sh
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=00:15:00
#SBATCH --job-name=load_counts
#SBATCH -o launch_scripts/cluster_logs/read_atac_counts_%a_%j.log
#SBATCH --array=1-500
#SBATCH --mem-per-cpu=1900MB

set -eo pipefail
module load R/4.0.5
TMPDIR=$(mktemp -d /tmp/$USER__$SLURM_JOBID-XXXX)

Rscript analysis/ATAC/MEGABULKS/0-extract_read_counts.R "$SLURM_ARRAY_TASK_ID"

echo "Tmp size is (after job run): $(du -sh $tempdir)"
rm -rf $tempdir

