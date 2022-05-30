#!/bin/sh
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=00:10:00
#SBATCH --job-name=analysis
#SBATCH -o launch_scripts/cluster_logs/snv_analysis.log 
#SBATCH --mem-per-cpu=5000MB

set -eo pipefail
module load R/4.0.5
TMPDIR=$(mktemp -d /tmp/$USER__$SLURM_JOBID-XXXX)

Rscript analysis/WGS/01-CNA/1-plot_cna.R

echo "Tmp size is (after job run): $(du -sh "$TMPDIR")"
rm -rf "$TMPDIR"
