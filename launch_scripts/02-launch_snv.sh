#!/bin/sh
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=80:00:00
#SBATCH --job-name=analysis
#SBATCH -o launch_scripts/cluster_logs/snv_analysis.log 
#SBATCH --mem-per-cpu=2500MB

set -eo pipefail
module load R/4.0.5
TMPDIR=$(mktemp -d /tmp/$USER__$SLURM_JOBID-XXXX)

#Start the job here.
sbatch --wait launch_scripts/02.1-load_vcfs.sh
Rscript analysis/WGS/02-SNV/1-load_snv_data.R
Rscript analysis/WGS/02-SNV/2-subset_driver_data.R
Rscript analysis/WGS/02-SNV/3-plot_drivers.R
Rscript analysis/WGS/02-SNV/4-plot_mutation_data.R 
sbatch --wait launch_scripts/02.2-lift_dnds.sh # only lifts vcfs to hg19
Rscript analysis/WGS/02-SNV/5.1-dnds.R
Rscript analysis/WGS/02-SNV/5.2-dnds.R
Rscript analysis/WGS/02-SNV/6-MSI_status.R

echo "Tmp size is (after job run): $(du -sh "$TMPDIR")"
rm -rf "$TMPDIR"
