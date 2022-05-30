#!/bin/sh
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=80:00:00
#SBATCH --job-name=tree_analysis
#SBATCH -o launch_scripts/cluster_logs/trees.log 
#SBATCH --mem-per-cpu=5000MB

set -eo pipefail
module load R/4.0.5
TMPDIR=$(mktemp -d /tmp/$USER__$SLURM_JOBID-XXXX)

# WGS trees:
Rscript analysis/WGS/03-TREES/1.1-reconstruct_mp_trees.R
Rscript analysis/WGS/03-TREES/1.2-plot_treesets.R
Rscript analysis/WGS/03-TREES/1.3-bootstrap_trees.R

# LP trees:
sbatch --wait ./launch_scripts/03.1-add_lp_samples.sh || Rscript analysis/WGS/03-TREES/2.1-assign_lowpass_samples_to_tree.R
Rscript analysis/WGS/03-TREES/2.2-plot_lowpass_trees.R
Rscript analysis/WGS/03-TREES/2.3-summary_lp_trees.R

# Annotated LP/WGS trees
analysis/WGS/03-TREES/3-plot_driver_annotated_tree.R

echo "Tmp size is (after job run): $(du -sh "$TMPDIR")"
rm -rf "$TMPDIR"
